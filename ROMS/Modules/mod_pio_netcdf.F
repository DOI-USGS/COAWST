#include "cppdefs.h"
      MODULE mod_pio_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
!
!git $Id$
!svn $Id: mod_pio_netcdf.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all Parallel-IO (PIO) variables definitions.   !
!  It also contains several variables and routines to facilitate the   !
!  manipulations of NetCDF data in parallel.  The PIO library was      !
!  developed at NCAR and it is the basis for the SCORPIO library.      !
!                                                                      !
!  For more information check:                                         !
!                                                                      !
!    https://e3sm.org/scorpio-parallel-io-library                      !
!                                                                      !
!  The libraries are available at GitHub:                              !
!                                                                      !
!    https://github.com/NCAR/ParallelIO                                !
!    https://github.com/E3SM-Project/scorpio                           !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE pio
!
      USE dateclock_mod, ONLY : datenum, datestr, time_units
      USE strings_mod,   ONLY : FoundError, find_string, lowercase
!
      implicit none
!
!  Interfaces for same name routine overloading. They differ in the kind
!  type parameter and data array rank.
!
      INTERFACE pio_netcdf_get_fatt              ! reads attributes
# ifdef SINGLE_PRECISION
        MODULE PROCEDURE pio_netcdf_get_fatt_dp
# endif
        MODULE PROCEDURE pio_netcdf_get_fatt_r8
      END INTERFACE pio_netcdf_get_fatt
!
      INTERFACE pio_netcdf_get_fvar              ! reads floating-point
# ifdef SINGLE_PRECISION
        MODULE PROCEDURE pio_netcdf_get_fvar_0dp
        MODULE PROCEDURE pio_netcdf_get_fvar_1dp
# endif
        MODULE PROCEDURE pio_netcdf_get_fvar_0d
        MODULE PROCEDURE pio_netcdf_get_fvar_1d
        MODULE PROCEDURE pio_netcdf_get_fvar_2d
        MODULE PROCEDURE pio_netcdf_get_fvar_3d
        MODULE PROCEDURE pio_netcdf_get_fvar_4d
      END INTERFACE pio_netcdf_get_fvar
!
      INTERFACE pio_netcdf_get_ivar              ! reads integer
        MODULE PROCEDURE pio_netcdf_get_ivar_0d
        MODULE PROCEDURE pio_netcdf_get_ivar_1d
        MODULE PROCEDURE pio_netcdf_get_ivar_2d
      END INTERFACE pio_netcdf_get_ivar
!
      INTERFACE pio_netcdf_get_lvar              ! reads logical
        MODULE PROCEDURE pio_netcdf_get_lvar_0d
        MODULE PROCEDURE pio_netcdf_get_lvar_1d
      END INTERFACE pio_netcdf_get_lvar

      INTERFACE pio_netcdf_get_svar              ! reads string
        MODULE PROCEDURE pio_netcdf_get_svar_0d
        MODULE PROCEDURE pio_netcdf_get_svar_1d
        MODULE PROCEDURE pio_netcdf_get_svar_2d
        MODULE PROCEDURE pio_netcdf_get_svar_3d
      END INTERFACE pio_netcdf_get_svar
!
      INTERFACE pio_netcdf_get_satt              ! gets string attribute
        MODULE PROCEDURE pio_netcdf_get_satt_g     ! global
        MODULE PROCEDURE pio_netcdf_get_satt_v     ! variable
      END INTERFACE pio_netcdf_get_satt
!
      INTERFACE pio_netcdf_get_time              ! reads time variable
        MODULE PROCEDURE pio_netcdf_get_time_0d
        MODULE PROCEDURE pio_netcdf_get_time_1d
      END INTERFACE pio_netcdf_get_time
!
      INTERFACE pio_netcdf_put_fvar              ! writes floating-point
# ifdef SINGLE_PRECISION
        MODULE PROCEDURE pio_netcdf_put_fvar_0dp
        MODULE PROCEDURE pio_netcdf_put_fvar_1dp
        MODULE PROCEDURE pio_netcdf_put_fvar_2dp
# endif
        MODULE PROCEDURE pio_netcdf_put_fvar_0d
        MODULE PROCEDURE pio_netcdf_put_fvar_1d
        MODULE PROCEDURE pio_netcdf_put_fvar_2d
        MODULE PROCEDURE pio_netcdf_put_fvar_3d
        MODULE PROCEDURE pio_netcdf_put_fvar_4d
      END INTERFACE pio_netcdf_put_fvar
!
      INTERFACE pio_netcdf_put_ivar              ! writes integer
        MODULE PROCEDURE pio_netcdf_put_ivar_0d
        MODULE PROCEDURE pio_netcdf_put_ivar_1d
        MODULE PROCEDURE pio_netcdf_put_ivar_2d
      END INTERFACE pio_netcdf_put_ivar
!
      INTERFACE pio_netcdf_put_lvar              ! writes logical
        MODULE PROCEDURE pio_netcdf_put_lvar_0d
        MODULE PROCEDURE pio_netcdf_put_lvar_1d
        MODULE PROCEDURE pio_netcdf_put_lvar_2d
      END INTERFACE pio_netcdf_put_lvar
!
      INTERFACE pio_netcdf_put_svar              ! writes string
        MODULE PROCEDURE pio_netcdf_put_svar_0d
        MODULE PROCEDURE pio_netcdf_put_svar_1d
        MODULE PROCEDURE pio_netcdf_put_svar_2d
        MODULE PROCEDURE pio_netcdf_put_svar_3d
      END INTERFACE pio_netcdf_put_svar
!
      PUBLIC :: pio_netcdf_check_dim    ! checks dimensions
      PUBLIC :: pio_netcdf_check_var    ! checks variables
      PUBLIC :: pio_netcdf_close        ! closes file
      PUBLIC :: pio_netcdf_create       ! creates file
      PUBLIC :: pio_netcdf_enddef       ! ends definition mode
      PUBLIC :: pio_netcdf_get_dim      ! reads dimension names/values
      PUBLIC :: pio_netcdf_inq_var      ! inquires variables
      PUBLIC :: pio_netcdf_inq_varid    ! inquires variable descriptor
      PUBLIC :: pio_netcdf_open         ! opens file
      PUBLIC :: pio_netcdf_redef        ! puts file in definition mode
      PUBLIC :: pio_netcdf_sync         ! synchronize file to disk
!
!  Switch indicating successful PIO initialization.
!
      logical :: LpioInitialized = .FALSE.

!
!  Set PIO internal level of debug information. The default value is 0,
!  allowed values 0-6.
!
      integer :: pio_debug = 0
!
!  Switch to debug creating, opening, and closing of NetCDF file
!  descriptor to monitor 'Too many open files' error. In the Unix
!  environment there is a limit to the number of open files. Use
!  'ulimit -a' or 'ulimit -S -n' to check.
!
#ifdef CHECK_OPEN_FILES
      logical :: Ldebug_pioFile = .TRUE.
#else
      logical :: Ldebug_pioFile = .FALSE.
#endif
!
!  Switch to create CDF-5 type pNetCDF files (PIO_64BIT_DATA) or
!  64-bit offset NetCDF-3 files (PIO_64BIT_OFFSET).  The CDF-5 are
!  not portable with third-part processing tools like Matlab.  They
!  are only included for benchmarking purposes. They are not
!  recommended for regular usage.
!
      logical :: TypeCDF5
!
!  Number of PIO components. In coupled systems, we can have
!  atmosphere, ocean, and other components.
!
      integer :: NpioComps = 1
!
!  PIO components indices.
!
      integer :: IpioROMS = 1                ! ROMS component
!
!  PIO dedicated computational and I/O ranks.
!
      character (len=1024) :: CcompRanks
      character (len=1024) :: CioRanks
!
!  PIO system descriptor handle (NpioComps,Ngrids). It describes the
!  system hardware arrangement: number of processes, which ones are
!  computational processes, and which ones are I/O processes.
!
      TYPE (IOSystem_desc_t), allocatable, target :: pioSystem(:,:)
!
!  IO descriptor structures for tile decomposition of single and double
!  precision data. They are used for the mapping between computational
!  and I/O processes. It describes how the data is distributed on the
!  computer.
!
      TYPE (IO_desc_t), pointer :: ioDesc_sp_p2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v2dvar(:) ! (i,j)
# if defined SSH_TIDES || defined UV_TIDES
      TYPE (IO_desc_t), pointer :: ioDesc_sp_rtides(:) ! (i,j,NTC)
# endif
# ifdef SOLVE3D
#  ifdef SEDIMENT
      TYPE (IO_desc_t), pointer :: ioDesc_sp_b3dvar(:) ! (i,j,Nbed)
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
      TYPE (IO_desc_t), pointer :: ioDesc_sp_l3dvar(:) ! (i,j,Nbands)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_l4dvar(:) ! (i,j,N,Nbands)
#  endif
      TYPE (IO_desc_t), pointer :: ioDesc_sp_p3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_w3dvar(:) ! (i,j,0:N)
# endif
# ifdef ADJUST_BOUNDARY
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r2dobc(:) ! (ij,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u2dobc(:) ! (ij,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v2dobc(:) ! (ij,4,Nbrec)
#  ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r3dobc(:) ! (ij,N,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u3dobc(:) ! (ij,N,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v3dobc(:) ! (ij,N,4,Nbrec)
#  endif
# endif
# ifdef ADJUST_STFLUX
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r2dfrc(:) ! (i,j,Nfrec)
# endif
# ifdef ADJUST_WSTRESS
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u2dfrc(:) ! (i,j,Nfrec)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v2dfrc(:) ! (i,j,Nfrec)
# endif
!
      TYPE (IO_desc_t), pointer :: ioDesc_dp_p2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u2dvar(:) ! (i,j)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v2dvar(:) ! (i,j)
# if defined SSH_TIDES || defined UV_TIDES
      TYPE (IO_desc_t), pointer :: ioDesc_dp_rtides(:) ! (i,j,NTC)
# endif
# ifdef SOLVE3D
#  ifdef SEDIMENT
      TYPE (IO_desc_t), pointer :: ioDesc_dp_b3dvar(:) ! (i,j,Nbed)
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
      TYPE (IO_desc_t), pointer :: ioDesc_dp_l3dvar(:) ! (i,j,Nbands)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_l4dvar(:) ! (i,j,N,Nbands)
#  endif
      TYPE (IO_desc_t), pointer :: ioDesc_dp_p3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v3dvar(:) ! (i,j,N)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_w3dvar(:) ! (i,j,0:N)
# endif
# ifdef ADJUST_BOUNDARY
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r2dobc(:) ! (ij,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u2dobc(:) ! (ij,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v2dobc(:) ! (ij,4,Nbrec)
#  ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r3dobc(:) ! (ij,N,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u3dobc(:) ! (ij,N,4,Nbrec)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v3dobc(:) ! (ij,N,4,Nbrec)
#  endif
# endif
# ifdef ADJUST_STFLUX
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r2dfrc(:) ! (i,j,Nfrec)
# endif
# ifdef ADJUST_WSTRESS
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u2dfrc(:) ! (i,j,Nfrec)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v2dfrc(:) ! (i,j,Nfrec)
# endif
!
!  IO descriptor structures for tile decomposition of single and double
!  precision special perfect restart data. They are used for the mapping
!  between computational and I/O processes.
!
      TYPE (IO_desc_t), pointer :: ioDesc_sp_rubar(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_rvbar(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_rzeta(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_ubar(:)   ! (i,j,3)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_vbar(:)   ! (i,j,3)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_zeta(:)   ! (i,j,3)
# ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_sp_ruvel(:)  ! (i,j,0:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_rvvel(:)  ! (i,j,0:N,2)
#  if defined GLS_MIXING || defined MY25_MIXING
      TYPE (IO_desc_t), pointer :: ioDesc_sp_tkevar(:) ! (i,j,0:N,2)
#  endif
      TYPE (IO_desc_t), pointer :: ioDesc_sp_trcvar(:) ! (i,j,1:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_uvel(:)   ! (i,j,1:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_vvel(:)   ! (i,j,1:N,2)
# endif
!
      TYPE (IO_desc_t), pointer :: ioDesc_dp_rubar(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_rvbar(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_rzeta(:)  ! (i,j,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_ubar(:)   ! (i,j,3)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_vbar(:)   ! (i,j,3)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_zeta(:)   ! (i,j,3)
# ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_dp_ruvel(:)  ! (i,j,0:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_rvvel(:)  ! (i,j,0:N,2)
#  if defined GLS_MIXING || defined MY25_MIXING
      TYPE (IO_desc_t), pointer :: ioDesc_dp_tkevar(:) ! (i,j,0:N,2)
#  endif
      TYPE (IO_desc_t), pointer :: ioDesc_dp_trcvar(:) ! (i,j,1:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_uvel(:)   ! (i,j,1:N,2)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_vvel(:)   ! (i,j,1:N,2)
# endif

# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
!
!  IO descriptor structures for tile decomposition of single and double
!  detiding harmonics data. They are used for the mapping between
!  computational and I/O processes.
!
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r2dhar(:) ! (i,j,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u2dhar(:) ! (i,j,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v2dhar(:) ! (i,j,0:2*NTC)
#  ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_sp_r3dhar(:) ! (i,j,N,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_u3dhar(:) ! (i,j,N,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_v3dhar(:) ! (i,j,N,0:2*NTC)
#  endif
!
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r2dhar(:) ! (i,j,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u2dhar(:) ! (i,j,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v2dhar(:) ! (i,j,0:2*NTC)
#  ifdef SOLVE3D
      TYPE (IO_desc_t), pointer :: ioDesc_dp_r3dhar(:) ! (i,j,N,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_u3dhar(:) ! (i,j,N,0:2*NTC)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_v3dhar(:) ! (i,j,N,0:2*NTC)
#  endif
# endif
#if defined PROPAGATOR && defined CHECKPOINTING
!
!  IO descriptor structures for tile decomposition of single and double
!  for the GST packed state data. They are used for the mapping between
!  computational and I/O processes.
!
      TYPE (IO_desc_t), pointer :: ioDesc_sp_Bvec(:)   ! (Mstate,NCV)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_resid(:)  ! (Mstate)
      TYPE (IO_desc_t), pointer :: ioDesc_sp_SworkD(:) ! (LworkD)
!
      TYPE (IO_desc_t), pointer :: ioDesc_dp_Bvec(:)   ! (Mstate,NCV)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_resid(:)  ! (Mstate)
      TYPE (IO_desc_t), pointer :: ioDesc_dp_SworkD(:) ! (LworkD)
#endif
!
!  Generic variable descriptors for current NetCDF file.
!
      TYPE (Var_desc_t),  pointer :: var_desc(:)
!
!  PIO supported methods of reading/writing NetCDF files
!
!    PIO_iotype_pnetcdf  => parallel read/write of PnetCDF (CDF-5)
!    PIO_iotype_netcdf   => serial   read/write of NetCDF3 (classic)
!    PIO_iotype_netcdf4c => parallel read/serial write of NetCDF4 (HDF5)
!    PIO_iotype_netcdf4p => parallel read/write of NETCDF4 (HDF5)
!
      integer pio_method
!
!  PIO supported method names:
!
!    PIO_iotype_pnetcdf  => 'PNETCDF'
!    PIO_iotype_netcdf   => 'NETCDF'
!    PIO_iotype_netcdf4c => 'NETCDF4C'   ! serial write data compression
!    PIO_iotype_netcdf4p => 'NETCDF4P'
!
      character (len=10) :: pio_MethodName
!
!  Number of coupled PIO components.
!
      integer :: NpioComp = 1
!
!  PIO number of processes used for I/O. If the I/O decomposition is
!  identical to the computational decomposition, "pio_NumIOtasks" is
!  equal to NtileI*NtileJ.  Typically, it is advantageous to define
!  the I/O decomposition in smaller number of processes for efficiency
!  and to avoid mpi-communication bottleneck.
!
      integer :: pio_NumIOtasks
!
!  PIO stride step in the mpi-rank between I/O tasks.
!
      integer :: pio_stride
!
!  PIO option that can be used to offset the first I/O task. The
!  default base is task 1.
!
      integer :: pio_base
!
!  In intra-communications mode, all processors in OCN_COMM_WORLD are
!  involved in computations. A subset or all processors does I/O (and
!  also computations). The "pio_NumIOtasks" and "pio_stride" variables
!  specify the total number of I/O tasks and the stride between them
!  with respect to the ROMS mpi-communicator object, OCN_COMM_WORLD.
!  The optional "pio_base" parameter is used to shift the first I/O
!  task away from the first computational task. This is often desirable
!  because the application first computational task usually has higher
!  memory requirements than other processes. If the MPI-processes are
!  spread over several hardware nodes, it is highly recommended to use
!  a value for PIO_STRIDE that scatters the I/O processes over all
!  nodes. Avoid all the I/O processes occupying the same node.
!
!  In the inter-communications (asynchronous) mode, the I/O tasks are
!  a disjointed set of dedicated I/O processes and do not  perform
!  computations. It is possible to have several groups of computational
!  units running separate models (coupling) where all the I/O data are
!  sent to dedicated processes. This interface is still experimental in
!  the SCORPIO library and not recommended for use at this time.
!
!  PIO number of mpi-aggregators to use in intra-communication mode
!  to improve mpi-collective I/O performance.
!
      integer :: pio_aggregator
!
!  PIO rearrangement method for moving data between computational and
!  I/O processes. It provides the ability to rearrange data between the
!  computational and I/O decompositions:
!
!    PIO_rearr_box    => Use a PIO internal box rearrangement
!    PIO_rearr_subset => Use a PIO internal sub-setting rearrangement
!
!  In the box method, data is rearranged from computational to I/O
!  processes in a continuous manner to the data ordering in the file.
!  Since the ordering of data between computational and I/O partitions
!  may be different, the rearrangement will require all-to-all MPI
!  communications. Also, notice that each computing tile may transfer
!  data to one or more I/O processes.
!
!  In the subset method, each I/O process is associated with a subset
!  of computing processes.  The computing tile sends its data to a
!  unique I/O process. The data on I/O processes may be more fragmented
!  to the ordering on disk, which may increase the communications to
!  the storage medium. However, the rearrangement scales better since
!  all-to-all MPI communications are not required.
!
      integer :: pio_rearranger
!
!  PIO rearranger communication between computational and I/O processes:
!
!    PIO_rearr_comm_p2p  => point-to-point communications
!    PIO_rearr_comm_coll => collective communications
!
      integer :: pio_rearr_comm
!
!  PIO rearranger communication betwen computational and I/O processes
!  flow control direction:
!
!    PIO_rearr_comm_fc_2d_disable => Disable flow control
!    PIO_rearr_comm_fc_2d_enable  => computational to IO processes, and
!                                    viceversa
!    PIO_rearr_comm_fc_1d_comp2io => computational to IO processes only
!    PIO_rearr_comm_fc_1d_io2comp => IO to computational processes only
!
!  Optimally, mpi-communications should be designed to send a modest
!  number messages evenly distributed accros a number of processes.  An
!  excessive number of messages to a single MPI-process can exhaust the
!  buffer space which can affect efficiency or failure due to the
!  slowdown in the retransmitting of dropped messages. PIO only send
!  messages (Isent) when the receiver is ready and has suffucient
!  resources.
!
      integer :: pio_rearr_fcd
!
!  PIO rearranger between computational to I/O processes (C2I):
!
      logical :: pio_rearr_C2I_HS                ! Enable handshake (T/F)
      logical :: pio_rearr_C2I_iS                ! Enable Isends (T/F)
!
      integer :: pio_rearr_C2I_PR                ! Max pending requests
!
!  PIO rearranger betwen I/O to computational processes (I2C):
!
      logical :: pio_rearr_I2C_HS                ! Enable handshake (T/F)
      logical :: pio_rearr_I2C_iS                ! Enable Isends (T/F)
!
      integer :: pio_rearr_I2C_PR                ! Max pending requests
!
!  External data representation for floating-point variables.
!
# ifdef OUT_DOUBLE
      integer, parameter :: PIO_FOUT = PIO_double
# else
      integer, parameter :: PIO_FOUT = PIO_real
# endif
# ifdef RST_SINGLE
      integer, parameter :: PIO_FRST = PIO_real
# else
      integer, parameter :: PIO_FRST = PIO_double
# endif
# ifdef DOUBLE_PRECISION
      integer, parameter :: PIO_TYPE = PIO_double
# else
      integer, parameter :: PIO_TYPE = PIO_real
# endif
!
!  External data representation for floating-point time and depth
!  variables. It is set to double precision for accuaracy in both
!  single and douple precision numerical kernel.
!
      integer, parameter :: PIO_TOUT = PIO_double
!
      CONTAINS
!
      SUBROUTINE pio_netcdf_copy_att (ng, model, VarName, AttName,      &
     &                                inp_ncname, inp_pioFile,          &
     &                                inp_VarID,                        &
     &                                out_ncname, out_pioFile,          &
     &                                out_VarID)
!
!=======================================================================
!                                                                      !
!  This function copy a global or variable attribute value from input  !
!  to output NetCDF file. It is done when output file is in define     !
!  mode.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     VarName      Variable name to process (string)                   !
!     AttName      Attribute name to process (string)                  !
!     inp_ncname   Input NetCDF filename (string)                      !
!     inp_pioFile  Input PIO file descriptor, TYPE(File_desc_t)        !
!                    inp_pioFile%fh         file handler               !
!                    inp_pioFile%iosystem   IO system descriptor       !
!     inp_VarID    Input NetCDF variable ID (integer)                  !
!     out_ncname   Output NetCDF filename (string)                     !
!     out_pioFile  Output PIO file descriptor, TYPE(File_desc_t)       !
!                    out_pioFile%fh         file handler               !
!                    out_pioFile%iosystem   IO system descriptor       !
!     out_VarID    Input NetCDF variable ID (integer)                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: inp_VarID, out_VarID
!
      character (len=*), intent(in) :: VarName, AttName
      character (len=*), intent(in) :: inp_ncname, out_ncname
!
      TYPE (File_desc_t), intent(inout) :: inp_pioFile, out_pioFile
!
!  Local variable declarations.
!
      logical :: foundit
!
      integer :: i, my_Atype, my_Natts
      integer :: status
      integer :: my_Aint(50)

      integer(PIO_OFFSET_KIND) :: my_Alen
!
      real(r8) :: my_Afloat(50)
!
      character (len=  40) :: my_Aname
      character (len=4096) :: my_Achar
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_copy_att"
!
      TYPE (Var_Desc_t) :: pioVar
!
!-----------------------------------------------------------------------
!  Copy global or variable attribute value from input to output NetCDF
!  file.
!-----------------------------------------------------------------------
!
!  Inquire about number of global or variable attributes.
!
      IF (inp_VarID.eq.PIO_global) THEN
        status=PIO_inquire(inp_pioFile,                                 &
     &                     nAttributes = my_Natts)
        IF (status.ne.PIO_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(inp_ncname),               &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        status=PIO_inq_varid(inp_pioFile, TRIM(VarName), pioVar)
        IF (status.eq.PIO_noerr) THEN
          status=PIO_inquire_variable(inp_pioFile, pioVar,              &
     &                                nAtts = my_Natts)
          IF (status.ne.PIO_noerr) THEN
            IF (Master) WRITE (stdout,20) TRIM(VarName),                &
     &                                    TRIM(inp_ncname),             &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(VarName),                  &
     &                                  TRIM(inp_ncname),               &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Inquire about requested global attribute value. If found, write its
!  value in the output NetCDF file.
!
      IF (status.eq.PIO_noerr) THEN
        DO i=1,my_Natts
          IF (inp_VarID.eq.PIO_global) THEN
            status=PIO_inq_attname(inp_pioFile, inp_VarID, i, my_Aname)
            IF (status.eq.PIO_noerr) THEN
              IF (TRIM(my_Aname).eq.TRIM(AttName)) THEN
                status=PIO_inq_att(inp_pioFile, inp_VarID,              &
     &                             TRIM(my_Aname),                      &
     &                             xtype = my_Atype,                    &
     &                             len = my_Alen)
                IF ((status.eq.PIO_noerr).and.                          &
     &              (my_Atype.eq.PIO_CHAR)) THEN
                  status=PIO_get_att(inp_pioFile, inp_VarID,            &
     &                               TRIM(my_Aname), my_Achar)
                  IF (status.eq.PIO_noerr) THEN
                    status=PIO_put_att(out_pioFile, PIO_global,         &
     &                                 TRIM(my_Aname), TRIM(my_Achar))
                    IF (FoundError(status, PIO_noerr,                   &
     &                  __LINE__, MyFile)) THEN
                      IF (Master) WRITE (stdout,40) TRIM(my_Aname),     &
     &                                              TRIM(out_ncname),   &
     &                                              TRIM(SourceFile)
                      exit_flag=3
                      ioerror=status
                      RETURN
                    END IF
                    EXIT
                  ELSE
                    IF (Master) WRITE (stdout,50) TRIM(my_Aname),       &
     &                                            TRIM(inp_ncname),     &
     &                                            TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                  END IF
                END IF
              END IF
            ELSE
              IF (Master) WRITE (stdout,60) 'global',                   &
     &                                      TRIM(inp_ncname),           &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
!
!  Inquire about requested variable attributes. Then, copy attributes
!  output NetCDF file variable.
!
          ELSE
            status=PIO_inq_attname(inp_pioFile, pioVar, i, my_Aname)
            IF (status.eq.PIO_noerr) THEN
              status=PIO_inq_att(inp_pioFile, pioVar,                   &
     &                           TRIM(my_Aname),                        &
     &                           xtype = my_Atype,                      &
     &                           len = my_Alen)
              IF (status.eq.PIO_noerr) THEN
                IF (my_Atype.eq.PIO_CHAR) THEN
                  status=PIO_get_att(inp_pioFile, pioVar,               &
     &                               TRIM(my_Aname),                    &
     &                               my_Achar)
                  IF (status.eq.PIO_noerr) THEN
                    status=PIO_put_att(out_pioFile, out_varID,          &
     &                                 TRIM(my_Aname),                  &
     &                                 TRIM(my_Achar))
                    IF (FoundError(status, PIO_noerr,                   &
     &                  __LINE__, MyFile)) THEN
                      IF (Master) WRITE (stdout,40) TRIM(my_Aname),     &
     &                                              TRIM(out_ncname),   &
     &                                              TRIM(SourceFile)
                      exit_flag=3
                      ioerror=status
                      RETURN
                    END IF
                  END IF
                ELSE IF (my_Atype.eq.PIO_INT) THEN
                  status=PIO_get_att(inp_pioFile, pioVar,               &
     &                               TRIM(my_Aname),                    &
     &                               my_Aint(1:my_Alen))
                  IF (status.eq.PIO_noerr) THEN
                    status=PIO_put_att(out_pioFile, out_VarID,          &
     &                                 TRIM(my_Aname),                  &
     &                                 my_Aint(1:my_Alen))
                    IF (FoundError(status, PIO_noerr,                   &
     &                  __LINE__, MyFile)) THEN
                      IF (Master) WRITE (stdout,40) TRIM(my_Aname),     &
     &                                              TRIM(out_ncname),   &
     &                                              TRIM(SourceFile)
                      exit_flag=3
                      ioerror=status
                      RETURN
                    END IF
                  END IF
                ELSE IF ((my_Atype.eq.PIO_REAL  ).or.                   &
     &                   (my_Atype.eq.PIO_DOUBLE)) THEN
                  status=PIO_get_att(inp_pioFile, pioVar,               &
     &                               TRIM(my_Aname),                    &
     &                               my_Afloat(1:my_Alen))
                  IF (status.eq.PIO_noerr) THEN
                    status=PIO_put_att(out_pioFile, out_VarID,          &
     &                                 TRIM(my_Aname),                  &
     &                                 my_Afloat(1:my_Alen))
                    IF (FoundError(status, PIO_noerr,                   &
     &                  __LINE__, MyFile)) THEN
                      IF (Master) WRITE (stdout,40) TRIM(my_Aname),     &
     &                                              TRIM(out_ncname),   &
     &                                              TRIM(SourceFile)
                      exit_flag=3
                      ioerror=status
                      RETURN
                    END IF
                  END IF
                END IF
              ELSE
                IF (Master) WRITE (stdout,50) TRIM(my_Aname),           &
     &                                        TRIM(inp_ncname),         &
     &                                        TRIM(SourceFile)
                exit_flag=2
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,60) 'variable',                 &
     &                                      TRIM(inp_ncname),           &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          END IF
        END DO
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while inquiring number ', &
     &        'of global attributes',/,23x,'in input file :',2x,a,      &
     &        /, 23x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while inquiring number ', &
     &        'attributes in variable :',2x,a,                          &
     &        /,23x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while inquiring ',        &
     &        'variable descriptor for:',2x,a,                          &
     &        /,23x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  40  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while writing ',          &
     &        'attribute:',2x,a,/,23x,'in output file:',2x,a,           &
     &        /,23x,'call from:',2x,a)
  50  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while reading ',          &
     &        'attribute:',2x,a,/,23x,'in input file:',2x,a,            &
     &        /,23x,'call from:',2x,a)
  60  FORMAT (/,' PIO_NETCDF_COPY_ATT - error while inquiring ',a,1x,   &
     &        'attributes',/,23x,'in input file:',2x,a,                 &
     &        /,23x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_copy_att
!
      FUNCTION pio_netcdf_find_var (ng, model, pioFile,                 &
     &                              VarName, pioVar) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This function checks if a requested variable is available in        !
!  a NetCDF file and returns its descriptor.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     VarName      Requested dimension name (string)                   !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     pioVar       PIO variable descriptor, TYPE(Var_Desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                    pioVar%vd          variable descriptor            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: VarName
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (Var_Desc_t),  intent(out)   :: pioVar
!
!  Local variable declarations.
!
      logical :: foundit
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_find_var"
!
!-----------------------------------------------------------------------
!  Inquire if requested variable is available in NetCDF file.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      status=PIO_noerr
      status=PIO_inq_varid(pioFile, TRIM(VarName), pioVar)
!
      IF (status.eq.PIO_noerr) THEN
        foundit=.TRUE.
      END IF
!
      RETURN
      END FUNCTION pio_netcdf_find_var
!
      SUBROUTINE pio_netcdf_get_dim (ng, model, ncname, pioFile,        &
     &                               DimName, DimSize, DimID)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition, if a particular dimension name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     DimName      Requested dimension name (string, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     DimSize      Size of requested dimension (integer, OPTIONAL)     !
!     DimID        ID of requested dimension (integer, OPTIONAL)       !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(out), optional :: DimSize
      integer, intent(out), optional :: DimID
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: DimName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical :: foundit
!
      integer :: i, j, status
      integer :: myID, myValue
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_dim"
!
      TYPE (File_desc_t) :: my_pioFile
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      ncformat=-1
      rec_id=-1
      rec_size=0
      dim_id=0
      dim_size=0
      DO i=1,Mdims
        DO j=1,LEN(dim_name(1))
          dim_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
         CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
         IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire file.
!
      status=PIO_inquire(my_pioFile, n_dim, n_var, n_gatt, rec_id)
      IF ((status.eq.PIO_noerr).and.(n_dim.le.Mdims)) THEN
!
!  Inquire about dimensions: names, ID, and size.
!
        rec_size=-1
        DO i=1,n_dim
          dim_id(i)=i
          status=PIO_inquire_dimension(my_pioFile, dim_id(i),           &
     &                                 dim_name(i), dim_size(i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            WRITE (stdout,10) dim_id(i), TRIM(ncname), TRIM(SourceFile)

            exit_flag=2
            ioerror=status
            EXIT
          END IF
          IF (dim_id(i).eq.rec_id) THEN
            rec_size=dim_size(i)
          END IF
        END DO
      ELSE
        IF (n_dim.gt.Mdims) THEN
          WRITE (stdout,20) ' Mdims = ', Mdims, n_dim
          exit_flag=2
          ioerror=0
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,30) TRIM(ncname), TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Load requested information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(DimName)) THEN
          DO i=1,n_dim
            IF (TRIM(dim_name(i)).eq.TRIM(DimName)) THEN
              foundit=.TRUE.
              myID=dim_id(i)
              myValue=dim_size(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(DimSize)) THEN
              DimSize=myValue
            END IF
            IF (PRESENT(DimID)) THEN
              DimID=myID
            END IF
          ELSE
            WRITE (stdout,40) TRIM(DimName), TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_DIM - error while reading dimension',  &
     &          ' ID:',2x,i0,/,22x,'in input file:',2x,a,               &
     &        /,22x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_DIM - too small dimension parameter,', &
     &        a,2i5,/,22x,'change file  mod_netcdf.F  and recompile')
  30  FORMAT (/,' PIO_NETCDF_GET_DIM - unable to inquire about',        &
     &          ' contents of input NetCDF file:',2x,a,                 &
     &        /,22x,'call from:',2x,a)
  40  FORMAT (/,' PIO_NETCDF_GET_DIM - requested dimension: ',a,/22x,   &
     &        'not found in input file:',2x,a,/,22x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_dim
!
      SUBROUTINE pio_netcdf_check_dim (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  It checks the file dimensions against model dimension parameters    !
!  for consitency.  All the dimensions information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO filename (string)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_check_dim"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_get_dim (ng, model, ncname)
      ELSE
        CALL pio_netcdf_get_dim (ng, model, ncname,                     &
     &                           pioFile = pioFile)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check dimensions for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_dim
        SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
          CASE ('xi_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_psi,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_psi,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_u,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_u,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_v,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_v,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('IorJ')
            IF (dim_size(i).ne.IOBOUNDS(ng)%IorJ) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%IorJ,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
# ifdef SOLVE3D
          CASE ('s_rho')
            IF (dim_size(i).ne.N(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), N(ng),         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
# endif
# ifndef RBL4DVAR_FCT_SENSITIVITY
#  ifdef ADJUST_BOUNDARY
          CASE ('obc_adjust')
            IF (dim_size(i).ne.Nbrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nbrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#  endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CASE ('frc_adjust')
            IF (dim_size(i).ne.Nfrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nfrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#  endif
# endif
        END SELECT
      END DO
!
 10   FORMAT (/,' PIO_NETCDF_CHECK_DIM - inconsistent size of',         &
     &          ' dimension: ',a,2x,2i5,/,24x,'in file: ',a)
!
      RETURN
      END SUBROUTINE pio_netcdf_check_dim
!
      SUBROUTINE pio_netcdf_check_var (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine inquires the NetCDF file variables and check if the    !
!  values of few of them are consitent with the parameters provided    !
!  in input scripts. All the variables information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: IDmod, Npts, i, ib, ic, j, j1, j2, status
      integer :: IvarS
!
      real(r8), parameter :: RoundOff = 1.0e-7_r8

      real(r8) :: FvarS, FvarV(50), VarVal
!
      character (len=80) :: text

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_check_var"
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF variables.
!-----------------------------------------------------------------------
!
!  Limit model identifier. The profiling is limited to iNLM, iTLM, iRPM,
!  and iADM.
!
      IF ((model.lt.1).or.(model.gt.4)) THEN
        IDmod=iNLM
      ELSE
        IDmod=model
      END IF
!
! Inquire about all variables.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_inq_var (ng, IDmod, ncname)
      ELSE
        CALL pio_netcdf_inq_var (ng, IDmod, ncname,                     &
     &                           pioFile = pioFile)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check several important variables for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(i))))

#ifdef SOLVE3D
          CASE ('Vtransform')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_ivar (ng, IDmod, ncname,              &
     &                                  var_name(i), IvarS)
            ELSE
              CALL pio_netcdf_get_ivar (ng, IDmod, ncname,              &
     &                                  var_name(i), IvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (IvarS.ne.Vtransform(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vtransform(ng),      &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Vstretching')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_ivar (ng, IDmod, ncname,              &
     &                                  var_name(i), IvarS)
            ELSE
              CALL pio_netcdf_get_ivar (ng, IDmod, ncname,              &
     &                                  var_name(i), IvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (IvarS.ne.Vstretching(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(var_name(i)),          &
     &                                      IvarS, Vstretching(ng),     &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('hc')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS)
            ELSE
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(hc(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, hc(ng),              &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_s')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS)
            ELSE
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(theta_s(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_s(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('theta_b')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS)
            ELSE
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

            IF (ABS(theta_b(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, theta_b(ng),         &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
          CASE ('Tcline')
            IF (.not.PRESENT(pioFile)) THEN
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS)
            ELSE
              CALL pio_netcdf_get_fvar (ng, IDmod, ncname,              &
     &                                  var_name(i), FvarS,             &
     &                                  pioFile = pioFile)
            END IF
            IF (FoundError(exit_flag, NoError,__LINE__, MyFile)) RETURN

            IF (ABS(Tcline(ng)-FvarS).gt.RoundOff) THEN
              IF (Master) WRITE (stdout,20) TRIM(var_name(i)),          &
     &                                      FvarS, Tcline(ng),          &
     &                                      TRIM(ncname)
              exit_flag=5
              EXIT
            END IF
#endif
#ifdef FOUR_DVAR
          CASE ('Hgamma')
            IF ((model.eq.5).or.(model.eq.10).or.(model.eq.11)) THEN
              IF (.not.find_string(var_name,n_var,'HgammaM',ic)         &
     &            .and.(model.eq.5).and.(NSA.eq.2)) THEN
                VarVal=Hgamma(2)
# ifdef ADJUST_BOUNDARY
              ELSE IF (.not.find_string(var_name,n_var,'HgammaB',ic)    &
     &                 .and.(model.eq.10)) THEN
                VarVal=Hgamma(3)
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
              ELSE IF (.not.find_string(var_name,n_var,'HgammaF',ic)    &
     &                 .and.(model.eq.11)) THEN
                VarVal=Hgamma(4)
# endif
              ELSE                   ! Backward compatability logic
                VarVal=Hgamma(1)     ! for a single Hgamma value
              END IF
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(VarVal-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, VarVal,            &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# ifdef WEAK_CONSTRAINT
          CASE ('HgammaM')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(2)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(2),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# ifdef ADJUST_BOUNDARY
          CASE ('HgammaB')
            IF (model.eq.10) THEN
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(3)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(3),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CASE ('HgammaF')
            IF (model.eq.11) THEN
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(Hgamma(4)-FvarS).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Hgamma(4),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
# endif
# ifdef SOLVE3D
          CASE ('Vgamma')
            IF ((model.eq.5).or.(model.eq.10)) THEN
              IF (.not.find_string(var_name,n_var,'VgammaM',ic)         &
     &            .and.(model.eq.5).and.(NSA.eq.2)) THEN
                VarVal=Vgamma(2)
# ifdef ADJUST_BOUNDARY
              ELSE IF (.not.find_string(var_name,n_var,'VgammaB',ic)    &
     &                 .and.(model.eq.10)) THEN
                VarVal=Vgamma(3)
# endif
              ELSE                   ! Backward compatability logic
                VarVal=Vgamma(1)     ! for a single Vgamma value
              END IF
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-VarVal).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, VarVal,            &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  ifdef WEAK_CONSTRAINT
          CASE ('VgammaM')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-Vgamma(2)).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Vgamma(2),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  endif
#  ifdef ADJUST_BOUNDARY
          CASE ('VgammaB')
            IF (model.eq.5) THEN
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS)
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarS,           &
     &                                    pioFile = pioFile)
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (ABS(FvarS-Vgamma(3)).gt.RoundOff) THEN
                IF (Master) WRITE (stdout,20) TRIM(var_name(i)),        &
     &                                        FvarS, Vgamma(3),         &
     &                                        TRIM(ncname)
                exit_flag=5
                EXIT
              END IF
            END IF
#  endif
# endif
          CASE ('Hdecay')
            IF ((model.eq.5).or.(model.eq.11)) THEN
              Npts=UBOUND(Hdecay,DIM=2)
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              END IF

              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (model.eq.5) THEN
                j1=1                   ! first state variable to check
# ifdef SOLVE3D
                j2=isTvar(NT(ng))      ! last  state variable to check
# else
                j2=3                   ! last  state variable to check
# endif
              ELSE IF (model.eq.11) THEN
# ifdef SOLVE3D
                j1=isTvar(NT(ng))+1    ! first state variable to check
# else
                j1=4                   ! first state variable to check
# endif
                j2=Npts                ! last  state variable to check
              END IF

              DO j=j1,j2
                IF (ABS(Hdecay(1,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(1,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Hdecay(1,j,ng),         &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
# ifdef SOLVE3D
          CASE ('Vdecay')
            IF ((model.eq.5).or.(model.eq.11)) THEN
              Npts=UBOUND(Vdecay,DIM=2)
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              IF (model.eq.5) THEN
                j1=1                   ! first state variable to check
                j2=isTvar(NT(ng))      ! last  state variable to check
              ELSE IF (model.eq.11) THEN
                j1=isTvar(NT(ng))+1    ! first state variable to check
                j2=Npts                ! last  state variable to check
              END IF

              DO j=j1,j2
                IF (ABS(Vdecay(1,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(1,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Vdecay(1,j,ng),         &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
# endif
# ifdef WEAK_CONSTRAINT
          CASE ('HdecayM')
            IF ((model.eq.5).and.(NSA.eq.2)) THEN
              Npts=UBOUND(Hdecay,DIM=2)
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              j1=1                     ! first state variable to check
#  ifdef SOLVE3D
              j2=isTvar(NT(ng))        ! last  state variable to check
#  else
              j2=3                     ! last  state variable to check
#  endif
              DO j=j1,j2
                IF (ABS(Hdecay(NSA,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(2,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Hdecay(NSA,j,ng),       &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
#  ifdef SOLVE3D
          CASE ('VdecayM')
            IF ((model.eq.5).and.(NSA.eq.2)) THEN
              Npts=UBOUND(Vdecay,DIM=2)
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1/),                &
     &                                    total = (/Npts/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              j1=1                     ! first state variable to check
              j2=isTvar(NT(ng))        ! last  state variable to check

              DO j=j1,j2
                IF (ABS(Vdecay(NSA,j,ng)-FvarV(j)).gt.RoundOff) THEN
                  text=TRIM(var_name(i))//                              &
     &                 '(2,'//TRIM(Vname(1,idSvar(j)))//')'
                  IF (Master) WRITE (stdout,20) TRIM(text),             &
     &                                          FvarV(j),               &
     &                                          Vdecay(NSA,j,ng),       &
     &                                          TRIM(ncname)
                  exit_flag=5
                  EXIT
                END IF
              END DO
            END IF
#  endif
# endif
# ifdef ADJUST_BOUNDARY
          CASE ('HdecayB')
            IF (model.eq.10) THEN
#  ifdef SOLVE3D
              Npts=isTvar(NT(ng))
#  else
              Npts=3
#  endif
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1,1/),              &
     &                                    total = (/Npts,4/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1,1/),              &
     &                                    total = (/Npts,4/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              ic=0
              j1=1                     ! first state variable to check
              j2=Npts                  ! last  state variable to check

              DO ib=1,4
                DO j=j1,j2
                  ic=ic+1
                  IF (Lobc(ib,j,ng)) THEN
                    WRITE (text,"(a,'(',i1,',',a,')')")                 &
     &                          TRIM(var_name(i)), ib,                  &
     &                          TRIM(Vname(1,idSvar(j)))
                    IF (ABS(HdecayB(j,ib,ng)-FvarV(ic)).gt.             &
     &                  RoundOff) THEN
                      IF (Master) WRITE (stdout,20) TRIM(text),         &
     &                                              FvarV(ic),          &
     &                                              HdecayB(j,ib,ng),   &
     &                                              TRIM(ncname)
                      exit_flag=5
                      EXIT
                    END IF
                  END IF
                END DO
              END DO
            END IF
#  ifdef SOLVE3D
          CASE ('VdecayB')
            IF (model.eq.10) THEN
              Npts=isTvar(NT(ng))
              IF (.not.PRESENT(pioFile)) THEN
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    start = (/1,1/),              &
     &                                    total = (/Npts,4/))
              ELSE
                CALL pio_netcdf_get_fvar (ng, IDmod, ncname,            &
     &                                    var_name(i), FvarV,           &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1,1/),              &
     &                                    total = (/Npts,4/))
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN

              ic=0
              j1=1                     ! first state variable to check
              j2=Npts                  ! last  state variable to check

              DO ib=1,4
                DO j=j1,j2
                  ic=ic+1
                  IF (Lobc(ib,j,ng)) THEN
                    WRITE (text,"(a,'(',i1,',',a,')')")                 &
     &                          TRIM(var_name(i)), ib,                  &
     &                          TRIM(Vname(1,idSvar(j)))
                    IF (ABS(VdecayB(j,ib,ng)-FvarV(ic)).gt.             &
     &                  RoundOff) THEN
                      IF (Master) WRITE (stdout,20) TRIM(text),         &
     &                                              FvarV(ic),          &
     &                                              VdecayB(j,ib,ng),   &
     &                                              TRIM(ncname)
                      exit_flag=5
                      EXIT
                    END IF
                  END IF
                END DO
              END DO
            END IF
#  endif
# endif
#endif
        END SELECT
      END DO

 10   FORMAT (/,' PIO_NETCDF_CHECK_VAR - inconsistent value of',        &
     &        ' variable: ',a,2x,2i5,/,24x,'in file: ',a)
 20   FORMAT (/,' PIO_NETCDF_CHECK_VAR - inconsistent value of',        &
     &        ' variable: ',a,2x,2(1pe14.6),/,24x,'in file: ',a)
!
      RETURN
      END SUBROUTINE pio_netcdf_check_var
!
      SUBROUTINE pio_netcdf_inq_var (ng, model, ncname, pioFile,        &
     &                               myVarName, SearchVar, pioVar,      &
     &                               nVarDim, nVarAtt)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition,  if a particular variable name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     myVarName    Requested variable name (string, OPTIONAL)          !
!     SearchVar    Switch used when searching a variable over          !
!                    multiple NetCDF files (logical, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!     nVarDim      Number of variable dimensions (integer, OPTIONAL)   !
!     nVarAtt      Number of variable attributes (integer, OPTIONAL)   !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!  If the OPTIONAL argument myVarName is provided, the following       !
!  information for requested variable is also stored:                  !
!                                                                      !
!     n_vatt       Number of variable attributes                       !
!     n_vdim       Number of variable dimensions                       !
!     var_kind     Variable external data type                         !
!     var_Aname    Variable attribute names (1:n_vatt)                 !
!     var_Achar    Variable string attribute values (1:n_vatt)         !
!     var_Afloat   Variable float attribute values (1:n_vatt)          !
!     var_Aint     Variable integer attribute values (1:n_vatt)        !
!     var_Dids     Variable dimensions ID (1:n_vdim)                   !
!     var_Dname    Variable dimensions name (1:n_vdim)                 !
!     var_Dsize    Variable dimensions size (1:n_vdim)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(out), optional :: SearchVar
!
      integer, intent(in) :: ng, model
      integer, intent(out), optional :: nVarDim
      integer, intent(out), optional :: nVarAtt
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: myVarName
!
      TYPE (File_desc_t), intent(in),  optional :: pioFile
      TYPE (Var_desc_t),  intent(out), optional :: pioVar
!
!  Local variable declarations.
!
      logical :: foundit, WriteError
!
      integer :: i, j, status
      integer :: att_id, my_Atype, my_id

      integer(pio_offset_kind) :: my_Alen
!
      real(r4) :: my_Afloat
      real(r8) :: my_Adouble
!
      character (len=1024) :: text

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_inq_var"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      att_kind=-1
      var_id=0
      var_natt=0
      var_flag=0
      var_type=0
      var_ndim=0
      var_dim=0
      status=PIO_noerr
      DO i=1,Matts
        DO j=1,LEN(att_name(1))
          att_name(i)(j:j)=' '
        END DO
      END DO
      DO i=1,Mvars
        DO j=1,LEN(var_name(1))
          var_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire NetCDF file.
!
      status=PIO_inquire(my_pioFile, n_dim, n_var, n_gatt, rec_id)
      IF ((status.eq.PIO_noerr).and.(n_var.le.Mvars)) THEN
!
!  Inquire global attribute names and their external data type.
!
        DO i=1,MIN(Matts,n_gatt)
          att_id=i
          status=PIO_inq_attname(my_pioFile, PIO_global, att_id,        &
     &                           att_name(i))
          IF (status.eq.PIO_noerr) THEN
            status=PIO_inq_att(my_pioFile, PIO_global,                  &
     &                          TRIM(att_name(i)),                      &
     &                          xtype = att_kind(i))
            IF (status.ne.PIO_noerr) THEN
              IF (Master) WRITE (stdout,10) i, TRIM(ncname),            &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          ELSE
            IF (Master) WRITE (stdout,10) i, TRIM(ncname),              &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
            EXIT
          END IF
        END DO
!
!  Inquire about variables: name, ID, dimensions, data type, and number
!  of attributes.
!
        IF (status.eq.PIO_noerr) THEN
          DO i=1,n_var
            var_id(i)=i
            var_flag(i)=1
            status=PIO_inquire_variable(my_pioFile, var_id(i),          &
     &                                  var_name(i), var_type(i),       &
     &                                  var_ndim(i), var_dim(:,i),      &
     &                                  var_natt(i))
            IF (status.eq.PIO_noerr) THEN
              status=PIO_inq_varid(my_pioFile, TRIM(var_name(i)),       &
     &                             var_desc(i))
              IF (status.eq.PIO_noerr) THEN
                DO j=1,MIN(NvarA,var_natt(i))
                  status=PIO_inq_attname(my_pioFile, var_desc(i), j,    &
     &                                   var_Aname(j))
                  IF (status.eq.PIO_noerr) THEN
                    IF (TRIM(var_Aname(j)).eq.'water_points'.and.       &
     &                  (var_ndim(i).gt.0)) THEN
                      var_flag(i)=-1
                    END IF
                  ELSE
                    IF (Master) WRITE (stdout,20) j, TRIM(var_name(i)), &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                END DO
              ELSE
                IF (Master) WRITE (stdout,30) TRIM(var_name(i)),        &
     &                                        TRIM(ncname),             &
     &                                        TRIM(SourceFile)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            ELSE
              IF (Master) WRITE (stdout,40) var_id(i), TRIM(ncname),    &
     &                                      TRIM(SourceFile)
              exit_flag=2
              ioerror=status
              EXIT
            END IF
          END DO
        END IF
      ELSE
        IF (n_var.gt.Mvars) THEN
          IF (Master) WRITE (stdout,50) 'Mvars = ', Mvars, n_var
          exit_flag=2
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,60) TRIM(ncname), TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
!  Load requested requested variable information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(myVarName)) THEN
          var_Dids=-1
          var_Dsize=0
          var_Aint=0
          var_Afloat=0.0_r8
          DO i=1,NvarA
            DO j=1,LEN(var_Aname(1))
              var_Aname(i)(j:j)=' '
            END DO
            DO j=1,LEN(var_Achar(1))
              var_Achar(i)(j:j)=' '
            END DO
          END DO
          DO i=1,NvarD
            DO j=1,LEN(var_Dname(1))
              var_Dname(i)(j:j)=' '
            END DO
          END DO
!
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(myVarName)) THEN
              foundit=.TRUE.
              my_id=var_id(i)
              n_vatt=var_natt(i)
              n_vdim=var_ndim(i)
              my_pioVar=var_desc(my_id)
              DO j=1,n_vdim
                var_Dids(j)=var_dim(j,i)
              END DO
              var_kind=var_type(i)
            END IF
          END DO
          IF (foundit) THEN
            IF (PRESENT(pioVar)) THEN
              pioVar=my_pioVar
            END IF
            IF (PRESENT(nVarDim)) THEN
              nVarDim=n_vdim
            END IF
            IF (PRESENT(nVarAtt)) THEN
              nVarAtt=n_vatt
            END IF
          END IF
!
!  If founded requested variable, inquire about is dimensions and
!  attributes.
!
          IF (foundit) THEN
            DO i=1,n_vdim
              status=PIO_inquire_dimension(my_pioFile, var_Dids(i),     &
     &                                     var_Dname(i), var_Dsize(i))
              IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
                IF (Master) WRITE (stdout,70) i, TRIM(myVarName),       &
     &                                        TRIM(ncname),             &
     &                                        TRIM(SourceFile)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
            IF (status.eq.PIO_noerr) THEN
              DO i=1,MIN(NvarA, n_vatt)
                status=PIO_inq_attname(my_pioFile, my_pioVar, i,        &
     &                                 var_Aname(i))
                IF (status.eq.PIO_noerr) THEN
                  status=PIO_inq_att(my_pioFile, my_pioVar,             &
     &                               TRIM(var_Aname(i)),                &
     &                               xtype = my_Atype,                  &
     &                               len = my_Alen)
                  IF (status.eq.PIO_noerr) THEN
                    IF ((my_Alen.eq.1).and.                             &
     &                  (my_Atype.eq.PIO_INT)) THEN
                      status=PIO_get_att(my_pioFile, my_pioVar,         &
     &                                   TRIM(var_Aname(i)),            &
     &                                   var_Aint(i))
                      IF (FoundError(status, PIO_noerr,                 &
     &                               __LINE__, MyFile)) THEN
                        IF (Master) WRITE (stdout,80) 'integer',        &
     &                                            TRIM(var_Aname(i)),   &
     &                                            TRIM(myVarName),      &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF ((my_Alen.eq.1).and.                        &
     &                       (my_Atype.eq.PIO_REAL)) THEN
                      status=PIO_get_att(my_pioFile, my_pioVar,         &
     &                                   TRIM(var_Aname(i)),            &
     &                                   my_Afloat)
                      IF (FoundError(status, PIO_noerr,                 &
     &                               __LINE__, MyFile)) THEN
                        IF (Master) WRITE (stdout,80) 'float',          &
     &                                            TRIM(var_Aname(i)),   &
     &                                            TRIM(myVarName),      &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
# ifdef SINGLE_PRECISION
                      var_Afloat(i)=my_Afloat
# else
                      var_Afloat(i)=REAL(my_Afloat, r8)
# endif

                    ELSE IF ((my_Alen.eq.1).and.                        &
     &                       (my_Atype.eq.PIO_DOUBLE)) THEN
                      status=PIO_get_att(my_pioFile, my_pioVar,         &
     &                                   TRIM(var_Aname(i)),            &
     &                                   my_Adouble)
                      IF (FoundError(status, PIO_noerr,                 &
     &                               __LINE__, MyFile)) THEN
                        IF (Master) WRITE (stdout,80) 'float',          &
     &                                            TRIM(var_Aname(i)),   &
     &                                            TRIM(myVarName),      &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
# ifdef SINGLE_PRECISION
                      var_Afloat(i)=REAL(my_Adouble, r4)
# else
                      var_Afloat(i)=my_Adouble
# endif

                    ELSE IF (my_Atype.eq.PIO_CHAR) THEN
                      status=PIO_get_att(my_pioFile, my_pioVar,         &
     &                                   TRIM(var_Aname(i)),            &
     &                                   text(1:my_Alen))
                      IF (FoundError(status, PIO_noerr,                 &
     &                               __LINE__, MyFile)) THEN
                        IF (Master) WRITE (stdout,80) 'string',         &
     &                                            TRIM(var_Aname(i)),   &
     &                                            TRIM(myVarName),      &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                      var_Achar(i)=text(1:my_Alen)
                    END IF
                  ELSE
                    IF (Master) WRITE (stdout,90) i, TRIM(myVarName),   &
     &                                            TRIM(ncname),         &
     &                                            TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                ELSE
                  IF (Master) WRITE (stdout,90) i, TRIM(myVarName),     &
     &                                          TRIM(ncname),           &
     &                                          TRIM(SourceFile)
                  exit_flag=4
                  ioerror=status
                  EXIT
                END IF
              END DO
            END IF
          END IF
!
!  Ignore error message if requested variable not found when searching
!  over multiple input NetCDF files.
!
          IF (PRESENT(SearchVar)) THEN
            SearchVar=foundit
            WriteError=.FALSE.
          ELSE
            WriteError=.TRUE.
          END IF
          IF (.not.foundit.and.WriteError) THEN
            IF (Master) WRITE (stdout,100) TRIM(myVarName),             &
     &                                     TRIM(ncname),                &
     &                                     TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while inquiring global',   &
     &        ' attribute: ',i2.2,/,22x,'in input file:',2x,a,          &
     &        /,22x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while inquiring attribute',&
     &        1x,i0,' for variable: ',a,/,22x,'in input file:',2x,a,    &
     &        /,22x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while inquiring variable', &
     &        ' descriptor for',2x,a,/,22x,'in input file:',2x,a,       &
     &        /,22x,'call from:',2x,a)
  40  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while inquiring variable', &
     &        ' ID:',2x,i0,/,22x,'in input file:',2x,a,                 &
     &        /,22x,'call from:',2x,a)
  50  FORMAT (/,' PIO_NETCDF_INQ_VAR - too small dimension parameter,', &
     &        1x,a,2i5,/,22x,'change file  mod_netcdf.F  and recompile')
  60  FORMAT (/,' PIO_NETCDF_INQ_VAR - unable to inquire about',        &
     &        ' contents of input NetCDF file:',2x,a,                   &
     &        /,22x,'call from:',2x,a)
  70  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while inquiring dimension',&
     &        1x,i0,' for variable:',2x,a,/,22x,'in input file:',2x,a,  &
     &        /,22x,'call from:',2x,a)
  80  FORMAT (/,' PIO_NETCDF_INQ_VAR - error while reading ',a,         &
     &        'attribute:',1x,a,' for variable ',a,/,22x,               &
     &        'in input file:',2x,a,/,22x,'call from:',2x,a)
  90  FORMAT (/,' PIO_NETCDF_INQ_VAR - unable to inquire name of ',     &
     &        'attribute ',i0,' for variable ',a,/,22x,                 &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a,/,18x,a)
 100  FORMAT (/,' PIO_NETCDF_INQ_VAR - requested variable:',2x,a,/22x,  &
     &        'not found in input file:',2x,a,/,22x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_inq_var
!
      SUBROUTINE pio_netcdf_inq_varid (ng, model, ncname, myVarName,    &
     &                                 pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine inquires the requested NetCDF variable descriptor.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Requested variable name (string)                    !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in)  :: pioFile
      TYPE (Var_desc_t),  intent(out) :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_inq_varid"
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
      status=PIO_inq_varid(pioFile, TRIM(myVarName), pioVar)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=3
        ioerror=status
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_INQ_VARID - error while inquiring ID',     &
     &        ' for variable:',2x,a,/,24x,'in input file:',2x,a,/,      &
     &        24x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_inq_varid

# ifdef SINGLE_PRECISION
!
      SUBROUTINE pio_netcdf_get_fatt_dp (ng, model, ncname, pioVar,     &
     &                                   AttName, AttValue, foundit,    &
     &                                   pioFile)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable double-precision attribute(s). !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!     AttName      Attribute name to read (string array)               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (double precision array)            !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      real(dp), intent(out) :: AttValue(:)
!
      TYPE (Var_desc_t),  intent(in) :: pioVar
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, natts, status
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fatt_dp"
!
      TYPE (File_desc_t) :: my_pioFile
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process and initialize.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=0.0_dp
      END DO
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire about requested attribute value.
!
      IF (pioVar%varID.eq.PIO_global) THEN
        status=PIO_inquire(my_pioFile,                                  &
     &                     nAttributes = my_natts)
      ELSE
        status=PIO_inquire_variable(my_pioFile, pioVar,                 &
     &                              name = my_Vname,                    &
     &                              nAtts = my_natts)
      END IF
      IF (status.eq.PIO_noerr) THEN
        DO j=1,my_natts
          status=PIO_inq_attname(my_pioFile, pioVar, j, my_Aname)
          IF (status.eq.PIO_noerr) THEN
            DO i=1,natts
              IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                status=PIO_get_att(my_pioFile, pioVar,                  &
     &                             TRIM(AttName(i)), AttValue(i))
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(AttName(i)),       &
     &                                          TRIM(my_Vname),         &
     &                                          TRIM(ncname),           &
     &                                          TRIM(SourceFile)
                  exit_flag=2
                  ioerror=status
                END IF
                foundit(i)=.TRUE.
                EXIT
              END IF
            END DO
          ELSE
            IF (Master) WRITE (stdout,20) j,                            &
     &                                    TRIM(my_Vname),               &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
            EXIT
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(my_Vname),                   &
     &                                TRIM(ncname),                     &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FATT_DP - error while reading ',       &
     &        'attribute:',1x,a,'for variable',1x,a,                    &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FATT_DP - error while inquiring ',     &
     &        'attribute:',1x,i2.2,'for variable',1x,a,                 &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_FATT_DP - error while inquiring ',     &
     &        'number of attributes for variable:',1x,a,                &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fatt_dp
# endif
!
      SUBROUTINE pio_netcdf_get_fatt_r8 (ng, model, ncname, pioVar,     &
     &                                   AttName, AttValue, foundit,    &
     &                                   pioFile)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable floating-point attribute(s).   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!     AttName      Attribute name to read (string array)               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (real array)                        !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      real(r8), intent(out) :: AttValue(:)
!
      TYPE (Var_desc_t),  intent(in) :: pioVar
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, natts, status
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fatt"
!
      TYPE (File_desc_t) :: my_pioFile
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process and initialize.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=0.0_r8
      END DO
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire about requested attribute value.
!
      IF (pioVar%varID.eq.PIO_global) THEN
        status=PIO_inquire(my_pioFile,                                  &
     &                     nAttributes = my_natts)
      ELSE
        status=PIO_inquire_variable(my_pioFile, pioVar,                 &
     &                              name = my_Vname,                    &
     &                              nAtts = my_natts)
      END IF
      IF (status.eq.PIO_noerr) THEN
        DO j=1,my_natts
          status=PIO_inq_attname(my_pioFile, pioVar, j, my_Aname)
          IF (status.eq.PIO_noerr) THEN
            DO i=1,natts
              IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                status=PIO_get_att(my_pioFile, pioVar,                  &
     &                             TRIM(AttName(i)), AttValue(i))
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(AttName(i)),       &
     &                                          TRIM(my_Vname),         &
     &                                          TRIM(ncname),           &
     &                                          TRIM(SourceFile)
                  exit_flag=2
                  ioerror=status
                END IF
                foundit(i)=.TRUE.
                EXIT
              END IF
            END DO
          ELSE
            IF (Master) WRITE (stdout,20) j,                            &
     &                                    TRIM(my_Vname),               &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
            EXIT
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(my_Vname),                   &
     &                                TRIM(ncname),                     &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FATT_R8 - error while reading ',       &
     &        'attribute:',1x,a,'for variable',1x,a,                    &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FATT_R8 - error while inquiring ',     &
     &        'attribute:',1x,i2.2,'for variable',1x,a,                 &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_FATT_R8 - error while inquiring ',     &
     &        'number of attributes for variable:',1x,a,                &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fatt_r8
!
      SUBROUTINE pio_netcdf_get_satt_g (ng, model, ncname, varid,       &
     &                                  AttName, AttValue, foundit,     &
     &                                  pioFile)
!
!=======================================================================
!                                                                      !
!  This routine gets requested global string attribute(s).             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     varid        Global attribute ID (integer, PIO_global)           !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (string array)                      !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: varid
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      character (len=*), intent(out) :: AttValue(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, natts, status
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_satt"
!
      TYPE (File_desc_t) :: my_pioFile
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process and initialize.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=' '
      END DO
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire about requested attribute value.
!
      status=PIO_inquire(my_pioFile,                                    &
     &                   nAttributes = my_natts)
      IF (status.eq.PIO_noerr) THEN
        DO j=1,my_natts
          status=PIO_inq_attname(my_pioFile, varid, j, my_Aname)
          IF (status.eq.PIO_noerr) THEN
            DO i=1,natts
              IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                status=PIO_get_att(my_pioFile, varid,                   &
     &                             TRIM(AttName(i)), AttValue(i))
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(AttName(i)),       &
     &                                          TRIM(my_Vname),         &
     &                                          TRIM(ncname),           &
     &                                          TRIM(SourceFile)
                  exit_flag=2
                  ioerror=status
                END IF
                foundit(i)=.TRUE.
                EXIT
              END IF
            END DO
          ELSE
            IF (Master) WRITE (stdout,20) j,                            &
     &                                    TRIM(my_Vname),               &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
            EXIT
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(my_Vname),                   &
     &                                TRIM(ncname),                     &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SATT_G - error while reading ',        &
     &        'attribute:',1x,a,'for variable',1x,a,                    &
     &        /,25x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_SATT_G - error while inquiring ',      &
     &        'attribute:',1x,i2.2,'for variable',1x,a,                 &
     &        /,25x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_SATT_G - error while inquiring',       &
     &        ' number of attributes for variable:',1x,a,               &
     &        /,25x,'in input file:',2x,a,/,19x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_satt_g
!
      SUBROUTINE pio_netcdf_get_satt_v (ng, model, ncname, pioVar,      &
     &                                  AttName, AttValue, foundit,     &
     &                                  pioFile)
!
!=======================================================================
!                                                                      !
!  This routine gets requested variable string attribute(s).           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t)           !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!     varid        Variable ID for variable attribute or               !
!                    NF90_GLOBAL for a global attribute (integer)      !
!     AttName      Attribute name to read (string array)               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     AttValue     Attribute value (string array)                      !
!     foundit      Switch (T/F) activated when the requested           !
!                    attribute is found (logical array)                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: AttName(:)
!
      logical, intent(out) :: foundit(:)
!
      character (len=*), intent(out) :: AttValue(:)
!
      TYPE (Var_desc_t),  intent(in) :: pioVar
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, j, my_natts, natts, status
!
      character (len=40) :: my_Aname
      character (len=40) :: my_Vname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_satt"
!
      TYPE (File_desc_t) :: my_pioFile
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
!  Get number of variable attributes to process and initialize.
!
      natts=UBOUND(AttName, DIM=1)
      DO i=1,natts
        foundit(i)=.FALSE.
        AttValue(i)=' '
      END DO
!
!  If appropriate, open file for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Inquire about requested attribute value.
!
      IF (pioVar%varID.eq.PIO_global) THEN
        status=PIO_inquire(my_pioFile,                                  &
     &                     nAttributes = my_natts)
      ELSE
        status=PIO_inquire_variable(my_pioFile, pioVar,                 &
     &                              name = my_Vname,                    &
     &                              nAtts = my_natts)
      END IF
      IF (status.eq.PIO_noerr) THEN
        DO j=1,my_natts
          status=PIO_inq_attname(my_pioFile, pioVar, j, my_Aname)
          IF (status.eq.PIO_noerr) THEN
            DO i=1,natts
              IF (TRIM(my_Aname).eq.TRIM(AttName(i))) THEN
                status=PIO_get_att(my_pioFile, pioVar,                  &
     &                             TRIM(AttName(i)), AttValue(i))
                IF (FoundError(status, PIO_noerr,                       &
     &                         __LINE__, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(AttName(i)),       &
     &                                          TRIM(my_Vname),         &
     &                                          TRIM(ncname),           &
     &                                          TRIM(SourceFile)
                  exit_flag=2
                  ioerror=status
                END IF
                foundit(i)=.TRUE.
                EXIT
              END IF
            END DO
          ELSE
            IF (Master) WRITE (stdout,20) j,                            &
     &                                    TRIM(my_Vname),               &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
            EXIT
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(my_Vname),                   &
     &                                TRIM(ncname),                     &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SATT_V - error while reading ',        &
     &        'attribute:',1x,a,'for variable',1x,a,                    &
     &        /,25x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_SATT_V - error while inquiring ',      &
     &        'attribute:',1x,i2.2,'for variable',1x,a,                 &
     &        /,25x,'in input file:',2x,a,/,23x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_SATT_V - error while inquiring',       &
     &        ' number of attributes for variable:',1x,a,               &
     &        /,25x,'in input file:',2x,a,/,19x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_satt_v
!
# ifdef SINGLE_PRECISION
!
      SUBROUTINE pio_netcdf_get_fvar_0dp (ng, model, ncname, myVarName, &
     &                                    A, pioFile, start, total,     &
     &                                    broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested double-precision scalar variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (double precision)             !
!     min_val      Read data minimum value (double precision, OPTIONAL)!
!     max_val      Read data maximum value (double precision, OPTIONAL)!
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(1))      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_0dp"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a double-precision scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, my_A)
          A=my_A(1)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_0DP - error while reading ',      &
     &        'variable:',2x,a,/,27x,'in input file:',2x,a,             &
     &        /,27x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_0DP - error while inquiring ',    &
     &        'descriptor for variable:',2x,a,/,27x,'in input file:',   &
     &        2x,a,/,27x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_0dp
!
      SUBROUTINE pio_netcdf_get_fvar_1dp (ng, model, ncname, myVarName, &
     &                                    A, pioFile, start, total,     &
     &                                    broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested double-precision 1D-array variable     !
!  from specified NetCDF file.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (double precision)           !
!     min_val      Read data minimum value (double precision, OPTIONAL)!
!     max_val      Read data maximum value (double precision, OPTIONAL)!
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(0:))     !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(:,1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
!
      integer :: i, status

      integer, dimension(1) :: Asize
!
      real(dp) :: Afactor, Aoffset, Aspval

      real(dp), parameter :: Aepsilon = 1.0E-8_r8

      real(dp), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_1dp"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a double-precision 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar,           &
     &                          AttName, AttValue, foundit,             &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_dp
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_dp
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO i=1,Asize(1)                       ! zero out missing values
          IF ((foundit(3).and.(ABS(A(i)-Aspval).lt.Aepsilon)).or.       &
     &        (.not.foundit(3).and.(ABS(A(i)).ge.ABS(Aspval)))) THEN
            A(i)=0.0_dp
          END IF
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_1DP - error while reading ',      &
     &        'variable:',2x,a,/,27x,'in input file:',2x,a,             &
     &        /,27x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_1DP - error while inquiring ',    &
     &        'descriptor for variable:',2x,a,/,27x,'in input file:',   &
     &        2x,a,/,27x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_1dp
# endif
!
      SUBROUTINE pio_netcdf_get_fvar_0d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total,      &
     &                                   broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point scalar variable from    !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(1))      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      real(r8), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, my_A)
          A=my_A(1)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_0D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_0d
!
      SUBROUTINE pio_netcdf_get_fvar_1d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total,      &
     &                                   broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 1D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(0:))     !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(:,1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
!
      integer :: i, status

      integer, dimension(1) :: Asize
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar,           &
     &                          AttName, AttValue, foundit,             &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO i=1,Asize(1)                       ! zero out missing values
          IF ((foundit(3).and.(ABS(A(i)-Aspval).lt.Aepsilon)).or.       &
     &        (.not.foundit(3).and.(ABS(A(i)).ge.ABS(Aspval)))) THEN
            A(i)=0.0_r8
          END IF
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_1D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_1d
!
      SUBROUTINE pio_netcdf_get_fvar_2d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total,      &
     &                                   broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 2D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(0:,:))   !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(0:,0:))  !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A(:,:,1))  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
!
      integer :: i, j, status

      integer, dimension(2) :: Asize
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
      END IF
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar,           &
     &                          AttName, AttValue, foundit,             &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO j=1,Asize(2)                       ! zero out missing values
          DO i=1,Asize(1)
            IF ((foundit(3).and.(ABS(A(i,j)-Aspval).lt.Aepsilon)).or.   &
     &          (.not.foundit(3).and.(ABS(A(i,j)).ge.ABS(Aspval)))) THEN
              A(i,j)=0.0_r8
            END IF
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=Afactor*A(i,j)
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              A(i,j)=A(i,j)+Aoffset
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_2D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_2d
!
      SUBROUTINE pio_netcdf_get_fvar_3d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total,      &
     &                                   broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 3D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL). It is    !
!                    ignored since PIO library broadcasts the values   !
!                    to all member in the group by default.            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
!
      integer :: i, j, k, status

      integer, dimension(3) :: Asize
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_3d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
      END IF
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar,           &
     &                          AttName, AttValue, foundit,             &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO k=1,Asize(3)                       ! zero out missing values
          DO j=1,Asize(2)
            DO i=1,Asize(1)
              IF ((foundit(3).and.                                      &
     &             (ABS(A(i,j,k)-Aspval).lt.Aepsilon)).or.              &
     &            (.not.foundit(3).and.                                 &
     &             (ABS(A(i,j,k)).ge.ABS(Aspval)))) THEN
                A(i,j,k)=0.0_r8
              END IF
            END DO
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=Afactor*A(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                A(i,j,k)=A(i,j,k)+Aoffset
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_3D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_3D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_3d
!
      SUBROUTINE pio_netcdf_get_fvar_4d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total,      &
     &                                   broadcast, min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 4D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!     broadcast    Switch to broadcast read values from root to all    !
!                    members of the communicator in distributed-       !
!                    memory applications (logical, OPTIONAL,           !
!                    default=TRUE)                                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 4D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_fvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: broadcast
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit
!
      integer :: i, j, k, l, status

      integer, dimension(4) :: Asize
!
      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), parameter :: Aepsilon = 1.0E-8_r8

      real(r8), dimension(3) :: AttValue
!
      character (len=12), dimension(3) :: AttName

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_fvar_4d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=total(1)
        Asize(2)=total(2)
        Asize(3)=total(3)
        Asize(4)=total(4)
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
        Asize(2)=UBOUND(A, DIM=2)
        Asize(3)=UBOUND(A, DIM=3)
        Asize(4)=UBOUND(A, DIM=4)
      END IF
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar,           &
     &                          AttName, AttValue, foundit,             &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=AttValue(1)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=AttValue(2)
        END IF

        IF (.not.foundit(3)) THEN
          Aspval=spval_check
        ELSE
          Aspval=AttValue(3)
        END IF

        DO l=1,Asize(4)                       ! zero out missing values
          DO k=1,Asize(3)
            DO j=1,Asize(2)
              DO i=1,Asize(1)
                IF ((foundit(3).and.                                    &
     &               (ABS(A(i,j,k,l)-Aspval).lt.Aepsilon)).or.          &
     &              (.not.foundit(3).and.                               &
     &               (ABS(A(i,j,k,l)).ge.ABS(Aspval)))) THEN
                  A(i,j,k,l)=0.0_r8
                END IF
              END DO
            END DO
          END DO
        END DO

        IF (foundit(1)) THEN                  ! scale data
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=Afactor*A(i,j,k,l)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO l=1,Asize(4)
            DO k=1,Asize(3)
              DO j=1,Asize(2)
                DO i=1,Asize(1)
                  A(i,j,k,l)=A(i,j,k,l)+Aoffset
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_FVAR_4D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_FVAR_4D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_fvar_4d
!
      SUBROUTINE pio_netcdf_get_lvar_0d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical scalar variable from specified !
!  NetCDF file.  The variable can be stored as an interger (0 or 1) or !
!  as a character ('T' or 'F').   Reading a character variable is very !
!  inefficient in parallel I/O.                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (logical)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_lvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_lvar (ng, iNLM, 'f.nc', 'VarName', A(1))      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      logical, intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: my_type, status

      integer :: AI
      integer, dimension(1) :: my_AI
!
      character (len=1) :: Achar(1)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_lvar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        status=PIO_inquire_variable(my_pioFile, my_pioVar,              &
     &                              xtype = my_type)
        IF (status.eq.PIO_noerr) THEN
          IF (my_type.eq.PIO_int) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=PIO_get_var(my_pioFile, my_pioVar, start, total,   &
     &                           my_AI)
              AI=my_AI(1)
            ELSE
              status=PIO_get_var(my_pioFile, my_pioVar, AI)
            END IF
            IF (status.eq.PIO_noerr) THEN
              IF (AI.eq.0) THEN
                A=.FALSE.
              ELSE
                A=.TRUE.
              END IF
            END IF
          ELSE IF (my_type.eq.PIO_char) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=PIO_get_var(my_pioFile, my_pioVar, start, total,   &
     &                           Achar)
            ELSE
              status=PIO_get_var(my_pioFile, my_pioVar, Achar)
            END IF
            IF (status.eq.PIO_noerr) THEN
              A=.FALSE.
              IF ((Achar(1).eq.'t').or.(Achar(1).eq.'T')) THEN
                A=.TRUE.
              END IF
            END IF
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL PIO_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_LVAR_0D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_LVAR_0D - error while inquiring ',     &
     &        'type for variable:',2x,a,/,26x,'in input file:',2x,a,    &
     &        /,26x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_LVAR_0D - error while inquiring ',     &
     &        ' descriptor for variable:',2x,a,                         &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE PIO_netcdf_get_lvar_0d
!
      SUBROUTINE pio_netcdf_get_lvar_1d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested logical 1D-array variable from         !
!  specified  NetCDF file.  The  variable  can be stored as an         !
!  interger (0 or 1) or  as a character ('T' or 'F').  Reading         !
!  a character variable is very inefficient in parallel I/O.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (logical)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_lvar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_lvar (ng, iNLM, 'f.nc', 'VarName', A(0:))     !
!    CALL pio_netcdf_get_lvar (ng, iNLM, 'f.nc', 'VarName', A(:,1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      logical, intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: i, my_type, status

      integer, dimension(SIZE(A,1)) :: AI
!
      character (len=1), dimension(SIZE(A,1)) :: Achar

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_lvar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        status=PIO_inquire_variable(my_pioFile, my_pioVar,              &
     &                              xtype = my_type)
        IF (status.eq.PIO_noerr) THEN
          IF (my_type.eq.PIO_int) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=PIO_get_var(my_pioFile, my_pioVar, start, total,   &
     &                           AI)
            ELSE
              status=PIO_get_var(my_pioFile, my_pioVar, AI)
            END IF
            IF (status.eq.PIO_noerr) THEN
              DO i=1,SIZE(A,1)
                IF (AI(i).eq.0) THEN
                  A(i)=.FALSE.
                ELSE
                  A(i)=.TRUE.
                END IF
              END DO
            END IF
          ELSE IF (my_type.eq.PIO_char) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=PIO_get_var(my_pioFile, my_pioVar, start, total,   &
     &                           Achar)
            ELSE
              status=PIO_get_var(my_pioFile, my_pioVar, Achar)
            END IF
            IF (status.eq.PIO_noerr) THEN
              DO i=1,SIZE(A,1)
                A(i)=.FALSE.
                IF ((Achar(i).eq.'t').or.(Achar(i).eq.'T')) THEN
                  A(i)=.TRUE.
                END IF
              END DO
            END IF
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,10) TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_LVAR_1D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_LVAR_1D - error while inquiring ',     &
     &        'type for variable:',2x,a,/,26x,'in input file:',2x,a,    &
     &        /,26x,'call from:',2x,a)
  30  FORMAT (/,' PIO_NETCDF_GET_LVAR_1D - error while inquiring ',     &
     &        ' descriptor for variable:',2x,a,                         &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_lvar_1d
!
      SUBROUTINE pio_netcdf_get_ivar_0d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer scalar variable from specified !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (integer)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(1))      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status

      integer, dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_ivar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, my_A)
          A=my_A(1)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_IVAR_0D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_IVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_ivar_0d
!
      SUBROUTINE pio_netcdf_get_ivar_1d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 1D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(0:))     !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(:,1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_ivar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in an integer 1D-array variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_IVAR_1D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_IVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,25x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_ivar_1d
!
      SUBROUTINE pio_netcdf_get_ivar_2d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 2D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A)         !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(0:,:))   !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(0:,0:))  !
!    CALL pio_netcdf_get_ivar (ng, iNLM, 'f.nc', 'VarName', A(:,:,1))  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      integer, intent(out) :: A(:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_ivar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in an integer 2D-array variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If file descriptor is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_IVAR_2D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_IVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,25x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_ivar_2d
!
      SUBROUTINE pio_netcdf_get_svar_0d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string scalar variable from specified  !
!  NetCDF file. The CDL of the scalar variable has one-dimension in    !
!  the NetCDF file for the number of characters:                       !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to read a scalar string use:                                        !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (string)                       !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=LEN(A)), dimension(1) :: my_A

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_svar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a string scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, my_A)
          A=my_A(1)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SVAR_0D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' pio_NETCDF_GET_SVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_svar_0d
!
      SUBROUTINE pio_netcdf_get_svar_1d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 1D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 1D-array         !
!  variable has two-dimensions in the NetCDF file, and the first       !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to read a single array element at location (i) use:                 !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_svar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a string 1D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SVAR_1D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' pio_NETCDF_GET_SVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_svar_1d
!
      SUBROUTINE pio_netcdf_get_svar_2d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 2D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 1D-array         !
!  variable has three-dimensions in the NetCDF file, and the first     !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2)                F90  !
!                                                                      !
!  to read a single array element at location (i,j) use:               !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_svar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a string 2D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SVAR_2D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' pio_NETCDF_GET_SVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_svar_2d
!
      SUBROUTINE pio_netcdf_get_svar_3d (ng, model, ncname, myVarName,  &
     &                                   A, pioFile, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 3D-array variable or array      !
!  element from specified NetCDF file. The CDL of the 1D-array         !
!  variable has four-dimensions in the NetCDF file, and the first      !
!  dimension is the number of characters:                              !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1, dim2, dim3)          F90  !
!                                                                      !
!  to read a single array element at location (i,j,k) use:             !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable or array element (string)    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:,:,:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_svar_3d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a string 3D-array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),    &
     &                                TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_SVAR_3D - error while reading ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' pio_NETCDF_GET_SVAR_3D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,                          &
     &        /,26x,'in input file:',2x,a,/,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_svar_3d
!
      SUBROUTINE pio_netcdf_get_time_0d (ng, model, ncname, myVarName,  &
     &                                   Rdate, A,                      &
     &                                   pioFile, start, total,         &
     &                                   min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested time scalar variable from specified    !
!  NetCDF file. If the "units" attribute of the form:                  !
!                                                                      !
!    'time-units since YYYY-MM-DD hh:mm:ss'                            !
!                                                                      !
!  is different than provided reference date "Rdate", it converts to   !
!  elapsed time since "Rdate".                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     Rdate        Reference date (real; [1] seconds, [2] days)        !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(1))    !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(in) :: Rdate(2)

      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical :: JulianOffset = .FALSE.
      logical :: Ldebug = .FALSE.

      logical, dimension(1) :: got_units
      logical, dimension(2) :: foundit
!
      integer :: ind, lstr, status
      integer :: year, month, day, hour, minutes
!
      real(dp) :: Afactor, Aoffset, my_Rdate(2), seconds
      real(dp) :: dnum_old, dnum_new

      real(dp), dimension(1) :: my_A
      real(r8), dimension(2) :: AttValue
!
      character (len=12) :: AttName(2)
      character (len=22) :: dstr_old, dstr_new
      character (len=40) :: UnitsAtt(1), UnitsValue(1)
      character (len=40) :: Units, Ustring

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_time_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, my_A)
          A=my_A(1)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                    TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar, AttName,  &
     &                          AttValue, foundit,                      &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=REAL(AttValue(1),dp)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=REAL(AttValue(2),dp)
        END IF

        IF (foundit(1)) THEN                  ! scale data
          A=Afactor*A
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          A=A+Aoffset
          IF (time_ref.eq.-2) JulianOffset=.TRUE.
        END IF
      END IF
!
!  Get time variable "units" attribute and convert to elapsed time
!  since reference date. If Julian Day Number (days or seconds) and
!  'add_offset' attribute,
!
      UnitsAtt(1)='units'

      CALL pio_netcdf_get_satt (ng, model, ncname, my_pioVar, UnitsAtt, &
     &                          UnitsValue, got_units,                  &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (got_units(1)) THEN
          Units=TRIM(lowercase(UnitsValue(1)))
          lstr=LEN_TRIM(Units)
          ind=INDEX(Units,'since')
          IF (ind.gt.0) THEN
            CALL time_units (TRIM(Units), year, month, day, hour,       &
     &                       minutes, seconds)
            CALL datenum (my_Rdate, year, month, day, hour, minutes,    &
     &                    seconds)
            IF (Rdate(1).ne.my_Rdate(1)) THEN
              Ustring=Units(1:ind-2)
              SELECT CASE (TRIM(Ustring))
                CASE ('second', 'seconds')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A
                    ELSE
                      dnum_old=my_Rdate(2)+A
                    END IF
                    CALL datestr (dnum_old, .FALSE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    A=A-Rdate(2)              ! 'add_offset' added above
                  ELSE
                    A=(my_Rdate(2)+A)-Rdate(2)
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(2)+A
                    CALL datestr (dnum_new, .FALSE., dstr_new)
                  END IF
                CASE ('day', 'days')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A
                    ELSE
                      dnum_old=my_Rdate(1)+A
                    END IF
                    CALL datestr (dnum_old, .TRUE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    A=A-Rdate(1)              ! 'add_offset' added above
                  ELSE
                    A=(my_Rdate(1)+A)-Rdate(1)
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(1)+A
                    CALL datestr (dnum_new, .TRUE., dstr_new)
                  END IF
              END SELECT
            END IF
          END IF
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_TIME_0D - error while reading',        &
     &        ' variable:',2x,a,/,26x,'in input file:',2x,a,            &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_TIME_0D - error while inquiring ID',   &
     &        ' for variable:',2x,a,/,26x,'in input file:',2x,a,        &
     &        /,22x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_time_0d
!
      SUBROUTINE pio_netcdf_get_time_1d (ng, model, ncname, myVarName,  &
     &                                   Rdate, A,                      &
     &                                   pioFile, start, total,         &
     &                                   min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested time 1D-array variable from specified  !
!  NetCDF file.  If the "units" attribute of the form:                 !
!                                                                      !
!    'time-units since YYYY-MM-DD hh:mm:ss'                            !
!                                                                      !
!  is different than provided reference date "Rdate", it converts to   !
!  elapsed time since "Rdate".                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    time variable name (string)                         !
!     Rdate        Reference date (real; [1] seconds, [2] days)        !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array time variable (real)                  !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:))   !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,1))  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      real(dp), intent(in) :: Rdate(2)

      real(dp), intent(out), optional :: min_val
      real(dp), intent(out), optional :: max_val

      real(dp), intent(out) :: A(:)
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
!
!  Local variable declarations.
!
      logical :: JulianOffset = .FALSE.
      logical :: Ldebug = .FALSE.

      logical, dimension(1) :: got_units
      logical, dimension(2) :: foundit
!
      integer :: i, ind, lstr, status
      integer :: year, month, day, hour, minutes

      integer, dimension(1) :: Asize
!
      real(dp) :: Afactor, Aoffset, my_Rdate(2), seconds
      real(dp) :: dnum_old, dnum_new

      real(r8), dimension(2) :: AttValue
!
      character (len=12) :: AttName(2)
      character (len=22) :: dstr_old, dstr_new
      character (len=40) :: UnitsAtt(1), UnitsValue(1)
      character (len=40) :: Units, Ustring

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_get_time_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a time 1D-array variable.
!-----------------------------------------------------------------------
!
      IF (PRESENT(start).and.PRESENT(total)) THEN
        Asize(1)=1
        DO i=1,SIZE(total)              ! this logic is for the case
          Asize(1)=Asize(1)*total(i)    ! of reading multidimensional
        END DO                          ! data into a compact 1D array
      ELSE
        Asize(1)=UBOUND(A, DIM=1)
      END IF
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  Read in time variable.
!
      status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
      IF (status.eq.PIO_noerr) THEN
        IF (PRESENT(start).and.PRESENT(total)) THEN
          status=PIO_get_var(my_pioFile, my_pioVar, start, total, A)
        ELSE
          status=PIO_get_var(my_pioFile, my_pioVar, A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),            &
     &                      TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      ELSE
        WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),              &
     &                    TRIM(SourceFile)
        exit_flag=2
        ioerror=status
      END IF
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '

      CALL pio_netcdf_get_fatt (ng, model, ncname, my_pioVar, AttName,  &
     &                          AttValue, foundit,                      &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (.not.foundit(1)) THEN
          Afactor=1.0_r8
        ELSE
          Afactor=REAL(AttValue(1),dp)
        END IF

        IF (.not.foundit(2)) THEN
          Aoffset=0.0_r8
        ELSE
          Aoffset=REAL(AttValue(2),dp)
        END IF

        IF (foundit(1)) THEN                  ! scale data
          DO i=1,Asize(1)
            A(i)=Afactor*A(i)
          END DO
        END IF

        IF (foundit(2)) THEN                  ! add data offset
          DO i=1,Asize(1)
            A(i)=A(i)+Aoffset
          END DO
          IF (time_ref.eq.-2) JulianOffset=.TRUE.
        END IF
      END IF
!
!  Get time variable "units" attribute and convert to elapsed time
!  since reference date.
!
      UnitsAtt(1)='units'

      CALL pio_netcdf_get_satt (ng, model, ncname, my_pioVar, UnitsAtt, &
     &                          UnitsValue, got_units,                  &
     &                          pioFile = my_pioFile)

      IF (exit_flag.eq.NoError) THEN
        IF (got_units(1)) THEN
          Units=TRIM(lowercase(UnitsValue(1)))
          lstr=LEN_TRIM(Units)
          ind=INDEX(Units,'since')
          IF (ind.gt.0) THEN
            CALL time_units (TRIM(Units), year, month, day, hour,       &
     &                       minutes, seconds)
            CALL datenum (my_Rdate, year, month, day, hour, minutes,    &
     &                    seconds)
            IF (Rdate(1).ne.my_Rdate(1)) THEN
              Ustring=Units(1:ind-2)
              SELECT CASE (TRIM(Ustring))
                CASE ('second', 'seconds')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A(1)
                    ELSE
                      dnum_old=my_Rdate(2)+A(1)
                    END IF
                    CALL datestr (dnum_old, .FALSE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    DO i=1,Asize(1)
                      A(i)=A(i)-Rdate(2)      ! 'add_offset' added above
                    END DO
                  ELSE
                    DO i=1,Asize(1)
                      A(i)=(my_Rdate(2)+A(i))-Rdate(2)
                    END DO
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(2)+A(1)
                    CALL datestr (dnum_new, .FALSE., dstr_new)
                  END IF
                CASE ('day', 'days')
                  IF (Ldebug) THEN
                    IF (JulianOffset) THEN
                      dnum_old=A(1)
                    ELSE
                      dnum_old=my_Rdate(1)+A(1)
                    END IF
                    CALL datestr (dnum_old, .TRUE., dstr_old)
                  END IF
                  IF (JulianOffset) THEN
                    DO i=1,Asize(1)
                      A(i)=A(i)-Rdate(1)      ! 'add_offset' added above
                    END DO
                  ELSE
                    DO i=1,Asize(1)
                      A(i)=(my_Rdate(1)+A(i))-Rdate(1)
                    END DO
                  END IF
                  IF (Ldebug) THEN
                    dnum_new=Rdate(1)+A(1)
                    CALL datestr (dnum_new, .TRUE., dstr_new)
                  END IF
              END SELECT
            END IF
          END IF
        END IF
      END IF
!
!  Compute minimum and maximum values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
!  If applicable, close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_GET_TIME_1D - error while reading',        &
     &        ' variable:',2x,a,/,26x,'in input file:',2x,a,            &
     &        /,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_GET_TIME_1D - error while inquiring ID',   &
     &        ' for variable:',2x,a,/,26x,'in input file:',2x,a,        &
     &        /,22x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_get_time_1d

#ifdef SINGLE_PRECISION
!
      SUBROUTINE pio_netcdf_put_fvar_0dp (ng, model, ncname, myVarName, &
     &                                    A, start, total,              &
     &                                    pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point scalar variable into a NetCDF  !
!  file.  If the file descritor is not provided, it opens the file,    !
!  writes data, and then closes the file.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(dp), intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_0dp"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
          status=PIO_put_var(my_pioFile, my_pioVar, A)
        ELSE
          my_A(1)=A
          status=PIO_put_var(my_pioFile, my_pioVar, start, total, my_A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_0DP - error while inquiring ',    &
     &        'descriptor for variable:',2x,a,/,27x,'in input file:',   &
     &        2x,a,/,27x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_0DP - error while writing ',      &
     &        'variable:',2x,a,/,27x,'in input file:',2x,a,             &
     &        /,27x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_0dp
!
      SUBROUTINE pio_netcdf_put_fvar_1dp (ng, model, ncname, myVarName, &
     &                                    A, start, total,              &
     &                                    pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 1D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(dp), intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_1dp"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_1DP - error while inquiring ',    &
     &        'descriptor for variable:',2x,a,/,27x,'in input file:',   &
     &        2x,a,/,27x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_1DP - error while writing ',      &
     &        'variable:',2x,a,/,27x,'in input file:',2x,a,             &
     &        /,27x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_1dp
!
      SUBROUTINE pio_netcdf_put_fvar_2dp (ng, model, ncname, myVarName, &
     &                                    A, start, total,              &
     &                                    pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 2D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (double precision)      !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(dp), intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_2dp"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_2DP - error while inquiring ',    &
     &        'descriptor for variable:',2x,a,/,27x,'in input file:',   &
     &        2x,a,/,27x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_2DP - error while writing ',      &
     &        'variable:',2x,a,/,27x,'in input file:',2x,a,             &
     &        /,27x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_2dp
#endif
!
      SUBROUTINE pio_netcdf_put_fvar_0d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point scalar variable into a NetCDF  !
!  file.  If the file descritor is not provided, it opens the file,    !
!  writes data, and then closes the file.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(r8), intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(r8), dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=pio_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
          status=PIO_put_var(my_pioFile, my_pioVar, A)
        ELSE
          my_A(1)=A
          status=PIO_put_var(my_pioFile, my_pioVar, start, total, my_A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_0D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_0d
!
      SUBROUTINE pio_netcdf_put_fvar_1d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 1D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(r8), intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_1D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_1d
!
      SUBROUTINE pio_netcdf_put_fvar_2d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 2D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(r8), intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_2D  - error while writing ',      &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_2d
!
      SUBROUTINE pio_netcdf_put_fvar_3d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 3D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(r8), intent(in) :: A(:,:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_3d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_3D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a,/,26x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_FVAR_3D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a,/,26x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_3d
!
      SUBROUTINE pio_netcdf_put_fvar_4d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a floating-point 4D-array variable into a NetCDF file.    !
!  If the file descritor is not provided, it opens the file, writes    !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      real(r8), intent(in) :: A(:,:,:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_fvar_4d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_FVAR_4D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a,/,26x,a)
  20  FORMAT (/,'PIO_NETCDF_PUT_FVAR_4D_PIO - error while writing ',    &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a,/,26x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_fvar_4d
!
      SUBROUTINE pio_netcdf_put_ivar_0d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It  writes an integer scalar variable into a NetCDF file. If the    !
!  file descritor is not provided, it opens the file, writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)

      integer, intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status

      integer, dimension(1) :: my_A
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_ivar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
          status=PIO_put_var(my_pioFile, my_pioVar, A)
        ELSE
          my_A(1)=A
          status=PIO_put_var(my_pioFile, my_pioVar, start, total, my_A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_IVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_IVAR_0D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_ivar_0d
!
      SUBROUTINE pio_netcdf_put_ivar_1d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes an integer 1D-array variable into a NetCDF file. If the   !
!  file descritor is not provided, it opens the file, writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)

      integer, intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_ivar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL PIO_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL PIO_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_IVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_IVAR_1D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_ivar_1d
!
      SUBROUTINE pio_netcdf_put_ivar_2d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes an integer 2D-array variable into a NetCDF file. If the   !
!  file descritor is not provided, it opens the file, writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)

      integer, intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pi_netcdf_put_ivar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_IVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_IVAR_2D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_ivar_2d
!
      SUBROUTINE pio_netcdf_put_lvar_0d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a logical scalar variable into a NetCDF file. If the file !
!  descritor is not provided, it opens the file, writes data, and then !
!  closes the file.                                                    !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO filename (string)                               !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      logical, intent(in) :: A
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
      integer :: AI

      integer, dimension(1) :: my_AI
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_lvar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
      IF (A) THEN
        AI=1
      ELSE
        AI=0
      END IF
!
!  Write out logical data as integers.
!
      IF (exit_flag.eq.NoError) THEN
        IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
          status=PIO_put_var(my_pioFile, my_pioVar, AI)
        ELSE
          my_AI(1)=AI
          status=PIO_put_var(my_pioFile, my_pioVar, start, total, my_AI)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_LVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,'PIO_NETCDF_PUT_LVAR_0D  - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_lvar_0d
!
      SUBROUTINE pio_netcdf_put_lvar_1d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a logical 1D-array variable into a NetCDF file. If the    !
!  file descritor is not provided, it opens the file, writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO filename (string)                               !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      logical, intent(in) :: A(:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: i, status

      integer, dimension(SIZE(A,1)) :: AI
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_lvar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
      DO i=1,SIZE(A,1)
        IF (A(i)) THEN
          AI(i)=1
        ELSE
          AI(i)=0
        END IF
      END DO
!
!  Write out logical data as integers.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, AI)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_LVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_LVAR_1D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_lvar_1d
!
      SUBROUTINE pio_netcdf_put_lvar_2d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  It writes a logical 2D-array variable into a NetCDF file. If the    !
!  file descritor is not provided, it opens the file, writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  The input logical data is converted to integer such that .FALSE.    !
!  is interpreted as zero, and .TRUE. is interpreted as one.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO filename (string)                               !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (logical)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: start(:), total(:)
!
      logical, intent(in) :: A(:,:)
!
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: i, j, status

      integer, dimension(SIZE(A,1),SIZE(A,2)) :: AI
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_lvar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
!  If file descriptor is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master)  WRITE (stdout,10) TRIM(myVarName), TRIM(ncname), &
     &                                   TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Convert logical data to integer: .FALSE. is interpreted as zero, and
!  .TRUE. is interpreted as one.
!
      DO j=1,SIZE(A,2)
        DO i=1,SIZE(A,1)
          IF (A(i,j)) THEN
            AI(i,j)=1
          ELSE
            AI(i,j)=0
          END IF
        END DO
      END DO
!
!  Write out logical data as integers.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, AI)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_LVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_LVAR_2D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_lvar_2d
!
      SUBROUTINE pio_netcdf_put_svar_0d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a string scalar variable into a file.  If       !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the scalar variable has        !
!  one-dimension in the NetCDF file for the number of characters:      !
!                                                                      !
!     char string(Nchars)                                         CDL  !
!                                                                      !
!     character (len=Nchars) :: string                            F90  !
!                                                                      !
!  to write a scalar string use:                                       !
!                                                                      !
!     start = (/1/)                                                    !
!     total = (/Nchars/)                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (string)                !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (1D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (1D vector integer)                     !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)
!
      character (len=*), intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=LEN(A)), dimension(1) :: my_A

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_svar_0d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Write out a scalar string.
!-----------------------------------------------------------------------
!
!  If file ID is not provided, open file for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        IF ((start(1).eq.1).and.(total(1).eq.1)) THEN
          status=PIO_put_var(my_pioFile, my_pioVar, A)
        ELSE
          my_A(1)=A
          status=PIO_put_var(my_pioFile, my_pioVar, start, total, my_A)
        END IF
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_SVAR_0D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_SVAR_0D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_svar_0d
!
      SUBROUTINE pio_netcdf_put_svar_1d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 1D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 1D-array variable has      !
!  two-dimensions in the NetCDF file, and the first dimension is       !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim1, Nchars)                                   CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1)                      F90  !
!                                                                      !
!  to write a single array element at location (i) use:                !
!                                                                      !
!     start = (/1, i/)                                                 !
!     total = (/Nchars, 1/)                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (1D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (2D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (2D vector integer)                     !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)
!
      character (len=*), intent(in) :: A(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_svar_1d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Write out a string 1D array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_SVAR_1D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_SVAR_1D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_svar_1d
!
      SUBROUTINE pio_netcdf_put_svar_2d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 2D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 1D-array variable has      !
!  three-dimensions in the NetCDF file, and the first dimension is     !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim2, dim1, Nchars)                             CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1,dim2)                 F90  !
!                                                                      !
!  to write a single array element at location (i,j) use:              !
!                                                                      !
!     start = (/1, i, j/)                                              !
!     total = (/Nchars, 1, 1/)                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (2D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (2D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (2D vector integer)                     !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)
!
      character (len=*), intent(in) :: A(:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_svar_2d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Write out a string 2D array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_SVAR_2D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_SVAR_2D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_svar_2d
!
      SUBROUTINE pio_netcdf_put_svar_3d (ng, model, ncname, myVarName,  &
     &                                   A, start, total,               &
     &                                   pioFile, pioVar)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 3D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file. The CDL of the 3D-array variable has      !
!  four-dimensions in the NetCDF file, and the first dimension is      !
!  the number of characters:                                           !
!                                                                      !
!     char string(dim3, dim2, dim1, Nchars)                       CDL  !
!                                                                      !
!     character (len=Nchars) :: string(dim1,dim2,dim3)            F90  !
!                                                                      !
!  to write a single array element at location (i,j,k) use:            !
!                                                                      !
!     start = (/1, i, j, k/)                                           !
!     total = (/Nchars, 1, 1, 1/)                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (3D string array)       !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (2D vector   !
!                    integer)                                          !
!     total        Number of data values to be written along each      !
!                    dimension (2D vector integer)                     !
!     pioFile      PIO file descriptor, TYPE(File_desc_t), OPTIONAL    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     pioVar       PIO variable descriptor, TYPE(Var_desc_t), OPTIONAL !
!                    pioVar%varID       Variable ID                    !
!                    pioVar%ncid        File ID                        !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)
!
      character (len=*), intent(in) :: A(:,:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
      TYPE (File_desc_t), intent(in), optional :: pioFile
      TYPE (Var_desc_t),  intent(in), optional :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_put_svar_3d"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Write out a string 3D array or array element.
!-----------------------------------------------------------------------
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_open (ng, model, TRIM(ncname), 1, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ELSE
        my_pioFile=pioFile
      END IF
!
!  If variable descriptor is not provided, inquire its value.
!
      IF (.not.PRESENT(pioVar)) THEN
        status=PIO_inq_varid(my_pioFile, TRIM(myVarName), my_pioVar)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      ELSE
        my_pioVar=pioVar
      END IF
!
!  Write out data.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_var(my_pioFile, my_pioVar, start, total, A)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),  &
     &                                  TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (.not.PRESENT(pioFile)) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_PUT_SVAR_3D - error while inquiring ',     &
     &        'descriptor for variable:',2x,a,/,26x,'in input file:',   &
     &        2x,a,/,26x,'call from:',2x,a)
  20  FORMAT (/,' PIO_NETCDF_PUT_SVAR_3D - error while writing ',       &
     &        'variable:',2x,a,/,26x,'in input file:',2x,a,             &
     &        /,26x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_put_svar_3d
!
      SUBROUTINE pio_netcdf_close (ng, model, pioFile, ncname, Lupdate)
!
!=======================================================================
!                                                                      !
!  This routine closes requested NetCDF file. If appropriate, it       !
!  also performs additional processing, like updating the global       !
!  attributes, before closing the file.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     ncname       NetCDF file name (string, OPTIONAL)                 !
!     Lupdate      Update global attribute (logical, OPTIONAl)         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      logical, intent(in), optional :: Lupdate
!
      character (len=*), intent(in), optional :: ncname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
!  Local variable declarations.
!
# ifdef BIOLOGY
      logical :: my_Lupdate
!
# endif
      integer :: FileH, i, status
# ifdef BIOLOGY
      integer :: is, ie, lstr
# endif
!
      character (len=200) :: my_ncname
# ifdef BIOLOGY
      character (len=512) :: bio_file
# endif

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_close"
!
!-----------------------------------------------------------------------
!  If open, close requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (ASSOCIATED(pioFile%iosystem)) THEN
        DO i=1,LEN(my_ncname)
          my_ncname(i:i)=' '
        END DO
!
        IF (.not.PRESENT(ncname)) THEN
!
!  Get filename, if any. It will be nice if there is a function in
!  the NetCDF library to do this. Fortunately, the filename is
!  written as a global attribute.
!
          status=PIO_get_att(pioFile, PIO_global, 'file', my_ncname)
        ELSE
          my_ncname=TRIM(ncname)
        END IF

# ifdef BIOLOGY
!
!  Determine updating value of biology header files global attribute.
!  This is only possible in output files. An error occurs in input
!  files open for reading only. This allows to use ROMS input files
!  with the "bio_file" attribute.
!
        IF (.not.PRESENT(Lupdate)) THEN
          my_Lupdate=.FALSE.
        ELSE
          my_Lupdate=Lupdate
        END IF
!
!  Update global attribute with the biology header files used.
!
        IF (my_Lupdate) THEN
          is=1
          DO i=1,512
            bio_file(i:i)=' '
          END DO
          DO i=1,4
            lstr=LEN_TRIM(BIONAME(i))
            IF (lstr.gt.0) THEN
              ie=is+lstr-1
              bio_file(is:ie)=TRIM(BIONAME(i))
              is=ie+1
              bio_file(is:is)=','
              is=is+2
            END IF
          END DO
          lstr=LEN_TRIM(bio_file)-1
          IF (lstr.gt.0) THEN
            status=PIO_put_att(pioFile, PIO_global, 'bio_file',         &
     &                         bio_file(1:lstr))
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) 'bio_file',                 &
     &                                      TRIM(my_ncname),            &
     &                                      TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
# endif
!
!  Close requested NetCDF file. After closing, the "pioFile%iosystem"
!  pointer becomes unassociated. Here, the "pioFile%fh" descriptor is
!  set to its closed state value of -1.
!
        IF (exit_flag.eq.NoError) THEN
          FileH=pioFile%fh
          IF (Ldebug_pioFile) THEN
            IF (Master) WRITE (DBout,'(a,1x," <= ",i8,2(2x,a))')        &
     &                        KernelString(model)//' PIO: CLOSE',       &
     &                        FileH, TRIM(my_ncname), TRIM(SourceFile)
            CALL my_flush (DBout)
          END IF
          CALL pio_closefile (pioFile)
!
          pioFile%fh=-1
        END IF
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_CLOSE - error while writing global ',      &
     &        'attribute:',2x,a,/,20x,'file:',2x,a,/,20x,               &
     &        'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_close
!
      SUBROUTINE pio_netcdf_create (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine creates a new NetCDF file. Currently, it only creates  !
!  file for serial or parallel I/O access.                             !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     Name of the new NetCDF file (string)                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pioVar     PIO variable descriptor, TYPE(Var_desc_t)             !
!                  pioVar%varID       Variable ID                      !
!                  pioVar%ncid        File ID                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(out) :: pioFile
!
!  Local variable declarations.
!
      integer :: iotype, my_cmode, status
!
      character (len=80) :: text

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_create"
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
!  The option PIO_64BIT_DATA create PnetCDF files of CDF-5 type. This
!  option is included for testing efficiency. However, CDF-5 files are
!  not portable in third party applications like Matlab. We recommend
!  Users to avoid creating CDF-5 files.
!
      SELECT CASE (pio_method)
        CASE (1)
          iotype=PIO_iotype_pnetcdf
          IF (TypeCDF5) THEN
            my_cmode=ior(PIO_CLOBBER, PIO_64BIT_DATA)
            text='Parallel read and write NetCDF3 file (CDF-5 type)'
          ELSE
            my_cmode=ior(PIO_CLOBBER, PIO_64BIT_OFFSET)
            text='Parallel read and write NetCDF3 file (64-bit offset)'
          END IF
        CASE (2)
          iotype=PIO_iotype_netcdf
          my_cmode=ior(PIO_CLOBBER, PIO_64BIT_OFFSET)
          text='Serial read and write NetCDF3 file (64-bit offset)'
        CASE (3)
          iotype=PIO_iotype_netcdf4c
          my_cmode=PIO_CLOBBER
          text='Parallel read and serial write NetCDF4/HDF5 file'
        CASE (4)
          iotype=PIO_iotype_netcdf4p
          my_cmode=PIO_CLOBBER
          text='Parallel read and write NetCDF4/HDF5 file'
      END SELECT
!
      status=PIO_CreateFile(pioSystem(IpioROMS,ng),                     &
     &                      pioFile,                                    &
     &                      iotype,                                     &
     &                      TRIM(ncname),                               &
     &                      my_cmode)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
        exit_flag=3
        ioerror=status
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(text), iotype
      END IF
!
      IF (Ldebug_pioFile) THEN
        IF (Master) WRITE (DBout,'(a," ** ",i8,2(2x,a))')               &
     &                    KernelString(model)//' PIO: CREATE',          &
     &                    pioFile%fh, TRIM(ncname), TRIM(SourceFile)
        CALL my_flush (DBout)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_CREATE - unable to create output NetCDF ', &
     &        'file:',/,21x,a,/,21x,'call from:',2x,a)
  20  FORMAT (21x,a,', ioType = ',i0)
!
      RETURN
      END SUBROUTINE pio_netcdf_create
!
      SUBROUTINE pio_netcdf_enddef (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine ends definition in an  open NetCDF  dataset.  The      !
!  changes made in define mode  are checked and committed to disk      !
!  if no errors occurred. The dataset is then placed in data mode,     !
!  so variable data can be read or written.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO file name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_enddef"
!
!-----------------------------------------------------------------------
!  Synchronize requested NetCDF file.
!-----------------------------------------------------------------------
!
      status=PIO_enddef(pioFile)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
        exit_flag=3
        ioerror=status
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_ENDDEF - unable to end definition mode',   &
     &        ' for file:',/,21x,a,/,21x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_enddef
!
      SUBROUTINE pio_netcdf_open (ng, model, ncname, omode, pioFILE)
!
!=======================================================================
!                                                                      !
!  This routine opens an existing NetCDF file for access. Currently,   !
!  it only opens file for serial or parallel I/O access.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     Name of the new NetCDF file (string)                  !
!     omode      Open mode flag:                                       !
!                    omode = 0,  read-only access (PIO_NOWRITE)        !
!                    omode = 1,  read and write access (PIO_WRITE)     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pioVar     PIO variable descriptor, TYPE(Var_desc_t)             !
!                  pioVar%varID       Variable ID                      !
!                  pioVar%ncid        File ID                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, omode
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(out) :: pioFile
!
!  Local variable declarations.
!
      integer :: my_omode, status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_open"
!
!-----------------------------------------------------------------------
!  Open requested NetCDF file.
!-----------------------------------------------------------------------
!
      SELECT CASE (omode)
        CASE (0)
          my_omode=PIO_nowrite
        CASE (1)
          my_omode=PIO_write
        CASE DEFAULT
          my_omode=PIO_nowrite
      END SELECT
!
      status=pio_openfile(pioSystem(IpioROMS,ng),                       &
     &                    pioFile,                                      &
     &                    pio_method,                                   &
     &                    ncname,                                       &
     &                    my_omode)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
        exit_flag=3
        ioerror=status
      END IF
!
      IF (Ldebug_pioFile) THEN
        IF (Master) WRITE (DBout,'(a,2x," => ",i8,2(2x,a))')            &
     &                    KernelString(model)//' PIO: OPEN',            &
     &                    pioFile%fh, TRIM(ncname), TRIM(SourceFile)
        CALL my_flush (DBout)
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_OPEN - unable to open existing NetCDF ',   &
     &        'file:',/,19x,a,/,19x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_open
!
      SUBROUTINE pio_netcdf_redef (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine puts an open NetCDF dataset into define mode, so       !
!  dimensions, variables, and attributes can be added or renamed       !
!  an attributes can be deleted.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
!  Local variable declarations.
!
      integer :: status
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_redef"
!
!-----------------------------------------------------------------------
!  Put open file into definition mode.
!-----------------------------------------------------------------------
!
      status=PIO_redef(pioFile)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
        exit_flag=3
        ioerror=status
      END IF
!
  10  FORMAT (/,' PIO_NETCDF_REDEF - unable to put in definition mode', &
     &        ' file:',/,20x,a,/,20x,'call from:',2x,a)
!
      RETURN
      END SUBROUTINE pio_netcdf_redef
!
      SUBROUTINE pio_netcdf_sync (ng, model, ncname, pioFile)
!
!=======================================================================
!                                                                      !
!  This routine synchronize to disk requested NetCDF file with         !
!  in-memory buffer to make data available to other  processes         !
!  immediately after it is written.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       PIO file name (string)                              !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t) :: pioFile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", pio_netcdf_sync"
!
!-----------------------------------------------------------------------
!  Synchronize requested NetCDF file.
!-----------------------------------------------------------------------
!
# if defined ASYNCHRONOUS_PIO || defined ASYNCHRONOUS_SCORPIO
      exit_flag=NoError
# else
      CALL PIO_syncfile (pioFile)
# endif
!
      RETURN
      END SUBROUTINE pio_netcdf_sync
#endif
      END MODULE mod_pio_netcdf
