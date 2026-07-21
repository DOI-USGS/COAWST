#include "cppdefs.h"
      MODULE set_pio_mod
#if defined PIO_LIB && defined DISTRIBUTE
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module contains several routines initializes and configures    !
!  Parallel-IO (PIO) in ROMS.                                          !
!                                                                      !
!=======================================================================
!
      USE pio
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
      USE mod_pio_netcdf
      USE mod_iounits
      USE mod_scalars
# if defined PROPAGATOR && defined CHECKPOINTING
      USE mod_storage
# endif
# ifdef VEGETATION
      USE mod_vegetation
# endif
!
# if defined SSH_TIDES || defined UV_TIDES
      USE mod_stepping, ONLY : NTC
# endif
      USE strings_mod,  ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: initialize_pio
      PUBLIC  :: finalize_pio
      PUBLIC  :: field_iodecomp
# if defined PROPAGATOR && defined CHECKPOINTING
      PUBLIC  :: state_iodecomp
# endif
      PUBLIC  :: set_iodecomp
!
      CONTAINS
!
      SUBROUTINE initialize_pio
!
!***********************************************************************
!                                                                      !
!  Initializes the PIO subsystem. It sets PIO decomposition for all    !
!  ROMS variables.                                                     !
!                                                                      !
!***********************************************************************
!
!  Local variable declarations.
!
# ifdef ASYNCHRONOUS_PIO
      logical, allocatable :: Lranks(:)
!
      integer :: ComputeSize

      integer, allocatable :: COMPUTE_COMM(:,:),  IO_COMM(:)
      integer, allocatable :: Compute_Ranks(:,:), IO_Ranks(:)
# endif
      integer :: MyError
      integer :: i, ic, ng
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", initialize_pio"
!
!-----------------------------------------------------------------------
!  Initialize PIO and get IO system descriptor.  It uses collective
!  communicatios.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(pioSystem)) THEN
        allocate ( pioSystem(NpioComps,Ngrids) )
      END IF
!
!  Set PIO internal level of debug information. The default value is 0,
!  allowed values 0-6.
!
      IF (pio_debug.gt.0) THEN
        CALL PIO_setdebuglevel (pio_debug)
      END IF

# if defined ASYNCHRONOUS_PIO
!
!-----------------------------------------------------------------------
!  If NCAR/UNIDATA PIO library, set the ranks of the computational and
!  IO processes with respect the initial (peer) communicator. Notice
!  that the communicator is natively split inside "PIO_init". Therefore,
!  the I/O processes do not return from the call. Instead go to an
!  internal loop and wait to receive further instructions from the
!  computational processes. I think that this is a better strategy.
!-----------------------------------------------------------------------
!
!  Set peer communicator as the initial ROMS communicator.
!
      PEER_COMM_WORLD=OCN_COMM_WORLD
      PeerSize=numthreads
      PeerRank=MyRank
      ComputeSize=PeerSize-pio_NumIOtasks
!
!  Allocate and initialize local arrays.
!
      IF (.not.allocated(COMPUTE_COMM)) THEN
        allocate ( COMPUTE_COMM(NpioComps,Ngrids) )
        COMPUTE_COMM=MPI_COMM_NULL
      END IF
      IF (.not.allocated(IO_COMM)) THEN
        allocate ( IO_COMM(Ngrids) )
        IO_COMM=MPI_COMM_NULL
      END IF
!
      IF (.not.allocated(Lranks)) THEN
        allocate ( Lranks(0:PeerSize-1) )
        Lranks(0:PeerSize-1)=.TRUE.
      END IF
      IF (.not.allocated(Compute_Ranks)) THEN
        allocate ( Compute_Ranks(ComputeSize,NpioComps) )
        Compute_Ranks=-1
      END IF
      IF (.not.allocated(IO_Ranks)) THEN
        allocate ( IO_Ranks(pio_NumIOtasks) )
        IO_Ranks=-1
      END IF
!
!  Set the ranks of the dedicated I/O processes with respect the peer
!  communicator.
!
      ic=pio_base
      IO_Ranks(1)=ic
      Lranks(ic)=.FALSE.
      DO i=2,pio_NumIOtasks
        IF ((ic+pio_stride).le.(PeerSize-1)) THEN
          ic=ic+pio_stride
          IO_Ranks(i)=ic
          Lranks(ic)=.FALSE.
        END IF
      END DO
      WRITE (CioRanks,'(*(i0,1x))') IO_Ranks
!
!  Set the ranks of the computational processes with respect the peer
!  communicator.
!
      ic=0
      DO i=0,PeerSize-1
        IF (Lranks(i)) THEN
          ic=ic+1
          Compute_Ranks(ic,IpioROMS)=i
        END IF
      END DO
      WRITE (CcompRanks,'(*(i0,1x))') Compute_Ranks
!
!  Create a new IO system for asynchronous or synchronous I/O.  The
!  asynchronous I/O is only possible PIO type ("io_pio") files.
!
!!    IF ((inp_lib.eq.io_pio).and.(out_lib.eq.io_pio)) THEN
      IF (out_lib.eq.io_pio) THEN
        DO ng=1,Ngrids
          CALL PIO_init (pioSystem(IpioROMS:,ng),                       &
     &                   PEER_COMM_WORLD,                               &
     &                   (/ComputeSize/),                               &
     &                   Compute_Ranks,                                 &
     &                   IO_Ranks,                                      &
     &                   pio_rearranger,                                &
     &                   COMPUTE_COMM(IpioROMS:,ng),                    &
     &                   IO_COMM(ng))
        END DO
!
!  Initialize ROMS kernel communicators.
!
        ng=1
        OCN_COMM_WORLD=COMPUTE_COMM(IpioROMS,ng)
        IO_COMM_WORLD =IO_COMM(ng)
!
!  Reset ROMS communicator parameters.
!
        CALL mpi_comm_size (OCN_COMM_WORLD, numthreads, MyError)
        CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
!
        Master=MyRank.eq.MyMaster
#  ifdef PARALLEL_IO
        InpThread=.TRUE.
        OutThread=.TRUE.
#  else
        IF (MyRank.eq.0) THEN
          InpThread=.TRUE.
          OutThread=.TRUE.
        ELSE
          InpThread=.FALSE.
          OutThread=.FALSE.
        END IF
#  endif
!
!  Deallocate local arrays.
!
        IF (allocated(COMPUTE_COMM))  deallocate (COMPUTE_COMM)
        IF (allocated(IO_COMM))       deallocate (IO_COMM)
        IF (allocated(Lranks))        deallocate (Lranks)
        IF (allocated(Compute_Ranks)) deallocate (Compute_Ranks)
        IF (allocated(IO_Ranks))      deallocate (IO_Ranks)
!
!  Otherwise, do synchronous I/O.
!
      ELSE
        DO ng=1,Ngrids
          CALL PIO_init (MyRank,                                        &
     &                   OCN_COMM_WORLD,                                &
     &                   pio_NumIOtasks,                                &
     &                   pio_aggregator,                                &
     &                   pio_stride,                                    &
     &                   pio_rearranger,                                &
     &                   pioSystem(IpioROMS,ng),                        &
     &                   base = pio_base)
        END DO
      END IF

# else
!
!-----------------------------------------------------------------------
!  Initialize synchronous PIO system.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL PIO_init (MyRank,                                          &
     &                 OCN_COMM_WORLD,                                  &
     &                 pio_NumIOtasks,                                  &
     &                 pio_aggregator,                                  &
     &                 pio_stride,                                      &
     &                 pio_rearranger,                                  &
     &                 pioSystem(IpioROMS,ng),                          &
     &                 base = pio_base)
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Set PIO rearrangement communication options.
!-----------------------------------------------------------------------
!
      LpioInitialized=.TRUE.
!
!  The rearranger communication type "pio_rearr_comm" has two choices:
!
!    PIO_rearr_comm_p2p               Point to point (send/recive)
!    PIO_rearr_comm_coll              Collective (gather/scatter)
!
!  The rearranger communication flow control direction "pio_rearr_fcd"
!  has four choices:
!
!    PIO_rearr_comm_fc_2d_enable      COMM to IO processes and viceversa
!    PIO_rearr_comm_fc_1d_comp2io     COMM to IO processes only
!    PIO_rearr_comm_fc_1d_io2comp     IO to COMM processes only
!    PIO_rearr_comm_fc_2d_disable     Disable flow control
!
!  Compute to IO (C2I) processes:
!
!    pio_rearr_C2I_HS                 Enable handshake (true/false)
!    pio_rearr_C2I_iS                 Enable Isends (true/false)
!    pio_rearr_C2I_PR                 Maximum pending requests
!
!  IO to compute (I2C) processes:
!
!    pio_rearr_I2C_HS                 Enable handshake (true/false)
!    pio_rearr_I2C_iS                 Enable Isends (true/false)
!    pio_rearr_I2C_PR                 Maximum pending requests
!
!  Use PIO_REARR_COMM_UNLIMITED_PEND_REQ for unlimited number of
!  requests.
!
      DO ng=1,Ngrids
        MyError=PIO_set_rearr_opts(pioSystem(IpioROMS,ng),              &
                                   pio_rearr_comm,                      &
     &                             pio_rearr_fcd,                       &
     &                             pio_rearr_C2I_HS,                    &
     &                             pio_rearr_C2I_iS,                    &
     &                             pio_rearr_C2I_PR,                    &
     &                             pio_rearr_I2C_HS,                    &
     &                             pio_rearr_I2C_iS,                    &
     &                             pio_rearr_I2C_PR)
        IF (FoundError(MyError, PIO_noerr, __LINE__, MyFile)) RETURN
      END DO
!
      RETURN
      END SUBROUTINE initialize_pio
!
      SUBROUTINE finalize_pio
!
!***********************************************************************
!                                                                      !
!  Finalizes the PIO subsystem.  It frees all the storage memory       !
!  associated with the IO decomposition.                               !
!                                                                      !
!***********************************************************************
!
!  Local variable declarations.
!
      integer :: i, ng, status
!
!-----------------------------------------------------------------------
!  Deallocate storage memory associated with IO decomposition.
!-----------------------------------------------------------------------
!
      IF (LpioInitialized) THEN
!
!  Single precision decomposition descriptors.
!
        DO ng=1,Ngrids
          DO i=1,NpioComps
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_p2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_r2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_u2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_v2dvar(ng))
# if defined SSH_TIDES || defined UV_TIDES
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_rtides(ng))
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_a3dvar(ng))
#  endif
#  ifdef SEDIMENT
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_b3dvar(ng))
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_l3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_l4dvar(ng))
#  endif
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_p3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_r3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_u3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_v3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_w3dvar(ng))

            IF (Nslice.gt.0) THEN
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_sp_p3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_sp_r3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_sp_u3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_sp_v3dvar(ng))
            END IF
# endif
          END DO
        END DO
!
!  Double precision decomposition descriptors.
!
        DO ng=1,Ngrids
          DO i=1,NpioComps
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_p2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_r2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_u2dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_v2dvar(ng))
# if defined SSH_TIDES || defined UV_TIDES
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_rtides(ng))
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_a3dvar(ng))
#  endif
#  ifdef SEDIMENT
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_b3dvar(ng))
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_l3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_l4dvar(ng))
#  endif
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_p3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_r3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_u3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_v3dvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_w3dvar(ng))

            IF (Nslice.gt.0) THEN
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_dp_p3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_dp_r3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_dp_u3dvar(ng))
             CALL PIO_freedecomp (pioSystem(i,ng), ioDesZ_dp_v3dvar(ng))
            END IF
# endif
          END DO
        END DO
!
!  Special restart and harmonics single precision decomposition
!  descriptors.
!
        DO ng=1,Ngrids
          DO i=1,NpioComps
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_rubar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_rvbar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_rzeta(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_ubar (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_vbar (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_zeta (ng))
# ifdef SOLVE3D
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_ruvel (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_rvvel (ng))
#  if defined GLS_MIXING || defined MY25_MIXING
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_tkevar(ng))
#  endif
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_trcvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_uvel  (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_vvel  (ng))
# endif
# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_r2dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_u2dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_v2dhar(ng))
#  ifdef SOLVE3D
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_r3dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_u3dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_v3dhar(ng))
#  endif
# endif
          END DO
        END DO
!
!  Special restart and harmonics double precision decomposition
!  descriptors.
!
        DO ng=1,Ngrids
          DO i=1,NpioComps
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_rubar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_rvbar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_rzeta(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_ubar (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_vbar (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_zeta (ng))
# ifdef SOLVE3D
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_ruvel (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_rvvel (ng))
#  if defined GLS_MIXING || defined MY25_MIXING
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_tkevar(ng))
#  endif
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_trcvar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_uvel  (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_vvel  (ng))
# endif
# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_r2dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_u2dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_v2dhar(ng))
#  ifdef SOLVE3D
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_r3dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_u3dhar(ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_v3dhar(ng))
#  endif
# endif
# if defined PROPAGATOR && defined CHECKPOINTING
!
!  I/O decomposition descriptors for GST single and double precision
!  state propagator data.
!
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_Bvec  (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_resid (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_sp_SworkD(ng))
!
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_Bvec  (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_resid (ng))
            CALL PIO_freedecomp (pioSystem(i,ng), ioDesc_dp_SworkD(ng))
# endif
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Shut down and clean up any memory associated with the PIO library.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO i=1,NpioComps
            CALL PIO_finalize (pioSystem(i,ng), status)
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE finalize_pio
!
      SUBROUTINE field_iodecomp (ng, ioSystem, ioType, ioDesc,          &
     &                           gtype, ndims, LBk, UBk, LBt, UBt)
!
!***********************************************************************
!                                                                      !
!  Sets the IO decomposition descriptor for ROMS field variable types. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     ioSystem     PIO system descriptor (TYPE IOSystem_desc_t)        !
!     ioType       PIO kind variable type (integer)                    !
!     gtype        Variable C-grid type (integer)                      !
!     ndims        Number of state variable dimensions (integer)       !
!     LBk          K- or 3rd-dimension Lower bound (integer, OPTIONAL) !
!     UBk          K- or 3rd-dimension Upper bound (integer, OPTIONAL) !
!     LBt          T- or 4th-dimension Lower bound (integer, OPTIONAL) !
!     UBt          T- or 4th-dimension Upper bound (integer, OPTIONAL) !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ioDesc       IO decomposition descriptor (TYPE io_desc_t)        !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ioType, gtype, ndims

      integer, intent(in), optional :: LBk, UBk
      integer, intent(in), optional :: LBt, UBt
!
      TYPE (IOSystem_desc_t), intent(in) :: ioSystem
      TYPE (io_desc_t), intent(out) :: ioDesc
!
!  Local variable declarations.
!
      integer :: Cgrid, ghost
      integer :: i, ic, j, jc, k, kc, l, lc, np
      integer :: Is, Ie, Js, Je
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ioff, Joff, Koff, Loff
      integer :: Ilen, Isize, Jlen, Jsize, Klen, Ksize, Llen, Lsize
      integer :: IJlen, IJKlen
      integer :: my_size
!
      integer(PIO_Offset_kind), allocatable :: map_decomp(:)
!
!-----------------------------------------------------------------------
!  Set the PIO computational decomposition for ROMS C-type variables
!  and array rank. It is based on variable kind type and its mapping
!  from storage order to memory order.
!-----------------------------------------------------------------------
!
!  Get GLOBAL lower and upper bounds for each variable type in input
!  or ouput NetCDF files.
!
      SELECT CASE (gtype)
        CASE (p2dvar, p3dvar)
          Cgrid=1
          Is=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Js=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, a3dvar, b3dvar, l3dvar, l4dvar, r3dvar)
          Cgrid=2
          Is=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Js=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Cgrid=3
          Is=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Js=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Cgrid=4
          Is=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Js=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE (w3dvar)
          Cgrid=2
          Is=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Js=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT
!
!  Get GLOBAL length for each variable dimension.
!
      Ilen=Ie-Is+1
      Jlen=Je-Js+1
      IJlen=Ilen*Jlen
!
      IF (PRESENT(LBk)) THEN
        IF (LBk.eq.0) THEN
          Koff=0
        ELSE
          Koff=1
        END IF
        Klen=UBk-LBk+1
        Ksize=Klen
        IJKlen=IJlen*Klen
      END IF
!
      IF (PRESENT(LBt)) THEN
        IF (LBt.eq.0) THEN
          Loff=0
        ELSE
          Loff=1
        END IF
        Llen=UBt-LBt+1
        Lsize=Llen
      END IF
!
!  Starting/ending I- and J-indices for each decomposition tile
!  according to C-grid locatation, excluding ghost points.
!
      ghost=0
      Imin=BOUNDS(ng) % Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng) % Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,MyRank)
!
!  Allocate 1D array for mapping of the storage order of the variable to
!  its memory order.
!
      Isize=Imax-Imin+1
      Jsize=Jmax-Jmin+1
!
      IF (ndims.eq.2) THEN
        my_size=Isize*Jsize
      ELSE IF (ndims.eq.3) THEN
        my_size=Isize*Jsize*Ksize
      ELSE IF (ndims.eq.4) THEN
        my_size=Isize*Jsize*Ksize*Lsize
      END IF
!
      IF (.not.ALLOCATED(map_decomp)) THEN
        allocate ( map_decomp(my_size) )
      END IF
      map_decomp=0_PIO_Offset_kind
!
!  Set variable decomposition mapping.
!
      IF (ndims.eq.2) THEN
        np=0
        DO j=Jmin,Jmax
          jc=(j-Joff)*Ilen
          DO i=Imin,Imax
            np=np+1
            ic=i+Ioff+jc
            map_decomp(np)=ic
          END DO
        END DO
      ELSE IF (ndims.eq.3) THEN
        np=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen+kc
            DO i=Imin,Imax
              np=np+1
              ic=i+Ioff+jc
              map_decomp(np)=ic
            END DO
          END DO
        END DO
      ELSE IF (ndims.eq.4) THEN
        np=0
        DO l=LBt,UBt
          lc=(l-Loff)*IJKlen
          DO k=LBk,UBk
            kc=(k-Koff)*IJlen+lc
            DO j=Jmin,Jmax
              jc=(j-Joff)*Ilen+kc
              DO i=Imin,Imax
                np=np+1
                ic=i+Ioff+jc
                map_decomp(np)=ic
              END DO
            END DO
          END DO
        END DO
      END IF
!
!  Set IO decomposition descriptor
!
      IF (ndims.eq.2) THEN
        CALL PIO_InitDecomp (ioSystem, ioType, (/Ilen,Jlen/),           &
     &                       map_decomp, ioDesc)
      ELSE IF (ndims.eq.3) THEN
        CALL PIO_InitDecomp (ioSystem, ioType, (/Ilen,Jlen,Klen/),      &
     &                       map_decomp, ioDesc)
      ELSE IF (ndims.eq.4) THEN
        CALL PIO_InitDecomp (ioSystem, ioType, (/Ilen,Jlen,Klen,Llen/), &
     &                       map_decomp, ioDesc)
      END IF
!
!  Deallocate.
!
      IF (allocated(map_decomp)) deallocate (map_decomp)
!
      RETURN
      END SUBROUTINE field_iodecomp
!
      SUBROUTINE set_iodecomp
!
!***********************************************************************
!                                                                      !
!  Sets the IO decomposition descriptors for ROMS input and output     !
!  variables.  They are used for the mapping between computational     !
!  and I/O processes.                                                  !
!                                                                      !
!***********************************************************************
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Allocate I/O decomposition descriptors.
!-----------------------------------------------------------------------
!
!  I/O decomposition descriptors for single precision data.
!
      allocate ( ioDesc_sp_p2dvar(Ngrids) )
# ifdef ADJUST_WSTRESS
      allocate ( ioDesc_sp_u2dfrc(Ngrids) )
      allocate ( ioDesc_sp_v2dfrc(Ngrids) )
# endif
      allocate ( ioDesc_sp_r2dvar(Ngrids) )
      allocate ( ioDesc_sp_u2dvar(Ngrids) )
      allocate ( ioDesc_sp_v2dvar(Ngrids) )
# if defined SSH_TIDES || defined UV_TIDES
      allocate ( ioDesc_sp_rtides(Ngrids) )
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
      allocate ( ioDesc_sp_a3dvar(Ngrids) )
#  endif
#  ifdef SEDIMENT
      allocate ( ioDesc_sp_b3dvar(Ngrids) )
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
      allocate ( ioDesc_sp_l3dvar(Ngrids) )
      allocate ( ioDesc_sp_l4dvar(Ngrids) )
#  endif
      allocate ( ioDesc_sp_p3dvar(Ngrids) )
#  ifdef ADJUST_STFLUX
      allocate ( ioDesc_sp_r2dfrc(Ngrids) )
#  endif
      allocate ( ioDesc_sp_r3dvar(Ngrids) )
      allocate ( ioDesc_sp_u3dvar(Ngrids) )
      allocate ( ioDesc_sp_v3dvar(Ngrids) )
      allocate ( ioDesc_sp_w3dvar(Ngrids) )

      IF (Nslice.gt.0) THEN                   ! output at depth slices
        allocate ( ioDesZ_sp_p3dvar(Ngrids) )
        allocate ( ioDesZ_sp_r3dvar(Ngrids) )
        allocate ( ioDesZ_sp_u3dvar(Ngrids) )
        allocate ( ioDesZ_sp_v3dvar(Ngrids) )
      END IF
# endif
!
! I/O decomposition descriptors for double precision data.
!
      allocate ( ioDesc_dp_p2dvar(Ngrids) )
# ifdef ADJUST_WSTRESS
      allocate ( ioDesc_dp_u2dfrc(Ngrids) )
      allocate ( ioDesc_dp_v2dfrc(Ngrids) )
# endif
      allocate ( ioDesc_dp_r2dvar(Ngrids) )
      allocate ( ioDesc_dp_u2dvar(Ngrids) )
      allocate ( ioDesc_dp_v2dvar(Ngrids) )
# if defined SSH_TIDES || defined UV_TIDES
      allocate ( ioDesc_dp_rtides(Ngrids) )
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
      allocate ( ioDesc_dp_a3dvar(Ngrids) )
#  endif
#  ifdef SEDIMENT
      allocate ( ioDesc_dp_b3dvar(Ngrids) )
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
      allocate ( ioDesc_dp_l3dvar(Ngrids) )
      allocate ( ioDesc_dp_l4dvar(Ngrids) )
#  endif
      allocate ( ioDesc_dp_p3dvar(Ngrids) )
#  ifdef ADJUST_STFLUX
      allocate ( ioDesc_dp_r2dfrc(Ngrids) )
#  endif
      allocate ( ioDesc_dp_r3dvar(Ngrids) )
      allocate ( ioDesc_dp_u3dvar(Ngrids) )
      allocate ( ioDesc_dp_v3dvar(Ngrids) )
      allocate ( ioDesc_dp_w3dvar(Ngrids) )

      IF (Nslice.gt.0) THEN                   ! output at depth slices
        allocate ( ioDesZ_dp_p3dvar(Ngrids) )
        allocate ( ioDesZ_dp_r3dvar(Ngrids) )
        allocate ( ioDesZ_dp_u3dvar(Ngrids) )
        allocate ( ioDesZ_dp_v3dvar(Ngrids) )
      END IF
# endif
!
!  I/O decomposition descriptors for special single precision
!  restart and harmonics data.
!
      allocate ( ioDesc_sp_rubar(Ngrids) )
      allocate ( ioDesc_sp_rvbar(Ngrids) )
      allocate ( ioDesc_sp_rzeta(Ngrids) )
      allocate ( ioDesc_sp_ubar (Ngrids) )
      allocate ( ioDesc_sp_vbar (Ngrids) )
      allocate ( ioDesc_sp_zeta (Ngrids) )
# ifdef SOLVE3D
      allocate ( ioDesc_sp_ruvel (Ngrids) )
      allocate ( ioDesc_sp_rvvel (Ngrids) )
#  if defined GLS_MIXING || defined MY25_MIXING
      allocate ( ioDesc_sp_tkevar(Ngrids) )
#  endif
      allocate ( ioDesc_sp_trcvar(Ngrids) )
      allocate ( ioDesc_sp_uvel  (Ngrids) )
      allocate ( ioDesc_sp_vvel  (Ngrids) )
# endif
# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
!
      allocate ( ioDesc_sp_r2dvar(Ngrids) )
      allocate ( ioDesc_sp_u2dvar(Ngrids) )
      allocate ( ioDesc_sp_v2dvar(Ngrids) )
#  ifdef SOLVE3D
      allocate ( ioDesc_sp_r3dvar(Ngrids) )
      allocate ( ioDesc_sp_u3dvar(Ngrids) )
      allocate ( ioDesc_sp_v3dvar(Ngrids) )
#  endif
# endif
!
!  I/O decomposition descriptors for special double precison
!  restart and harmonics data.
!
      allocate ( ioDesc_dp_rubar(Ngrids) )
      allocate ( ioDesc_dp_rvbar(Ngrids) )
      allocate ( ioDesc_dp_rzeta(Ngrids) )
      allocate ( ioDesc_dp_ubar (Ngrids) )
      allocate ( ioDesc_dp_vbar (Ngrids) )
      allocate ( ioDesc_dp_zeta (Ngrids) )
# ifdef SOLVE3D
      allocate ( ioDesc_dp_ruvel (Ngrids) )
      allocate ( ioDesc_dp_rvvel (Ngrids) )
#  if defined GLS_MIXING || defined MY25_MIXING
      allocate ( ioDesc_dp_tkevar(Ngrids) )
#  endif
      allocate ( ioDesc_dp_trcvar(Ngrids) )
      allocate ( ioDesc_dp_uvel  (Ngrids) )
      allocate ( ioDesc_dp_vvel  (Ngrids) )
# endif
# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
!
      allocate ( ioDesc_dp_r2dvar(Ngrids) )
      allocate ( ioDesc_dp_u2dvar(Ngrids) )
      allocate ( ioDesc_dp_v2dvar(Ngrids) )
#  ifdef SOLVE3D
      allocate ( ioDesc_dp_r3dvar(Ngrids) )
      allocate ( ioDesc_dp_u3dvar(Ngrids) )
      allocate ( ioDesc_dp_v3dvar(Ngrids) )
#  endif
# endif
# if defined PROPAGATOR && defined CHECKPOINTING
!
!  I/O decomposition descriptors for GST single and double precision
!  state propagator data. Its values are set in routine "wpoints_tile".
!
      allocate ( ioDesc_sp_Bvec(Ngrids) )
      allocate ( ioDesc_sp_resid(Ngrids) )
      allocate ( ioDesc_sp_SworkD(Ngrids) )
!
      allocate ( ioDesc_dp_Bvec(Ngrids) )
      allocate ( ioDesc_dp_resid(Ngrids) )
      allocate ( ioDesc_dp_SworkD(Ngrids) )
# endif
!
!-----------------------------------------------------------------------
!  Set the PIO computational decomposition for ROMS C-type variables
!  and array rank. It is based on variable kind type and its mapping
!  from storage order to memory order.
!-----------------------------------------------------------------------
!
!  Set I/O decomposition descriptors for single precision data
!
      DO ng=1,Ngrids
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_p2dvar(ng),                      &
     &                       p2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_r2dvar(ng),                      &
     &                       r2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_u2dvar(ng),                      &
     &                       u2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_v2dvar(ng),                      &
     &                       v2dvar, 2)
# if defined SSH_TIDES || defined UV_TIDES
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_rtides(ng),                      &
     &                       r2dvar+4, 3, 1, NTC(ng))
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_a3dvar(ng),                      &
     &                       a3dvar, 3, 1, NVEG)
#  endif
#  ifdef SEDIMENT
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_b3dvar(ng),                      &
     &                       b3dvar, 3, 1, Nbed)
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_l3dvar(ng),                      &
     &                       l3dvar, 3, 1, NDbands)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_l4dvar(ng),                      &
     &                       l4dvar, 4, 1, N(ng), NDbands)
#  endif
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_p3dvar(ng),                      &
     &                       p3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_r3dvar(ng),                      &
     &                       r3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_u3dvar(ng),                      &
     &                       u3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_v3dvar(ng),                      &
     &                       v3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_w3dvar(ng),                      &
     &                       w3dvar, 3, 0, N(ng))
!
        IF (Nslice.gt.0) THEN
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,    &
     &                         ioDesZ_sp_p3dvar(ng),                    &
     &                         p3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,    &
     &                         ioDesZ_sp_r3dvar(ng),                    &
     &                         r3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,    &
     &                         ioDesZ_sp_u3dvar(ng),                    &
     &                         u3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,    &
     &                         ioDesZ_sp_v3dvar(ng),                    &
     &                         v3dvar, 3, 1, Nslice)
        END IF
# endif
# if defined ADJUST_STFLUX && defined DISTRIBUTE
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_r2dfrc(ng),                      &
     &                       r2dvar, 3, 1, Nfrec(ng))
# endif
# ifdef ADJUST_WSTRESS
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_u2dfrc(ng),                      &
     &                       u2dvar, 3, 1, Nfrec(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_v2dfrc(ng),                      &
     &                       v2dvar, 3, 1, Nfrec(ng))
# endif
      END DO
!
!  Set IO decomposition descriptors for double precision data.
!
      DO ng=1,Ngrids
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_p2dvar(ng),                      &
     &                       p2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_r2dvar(ng),                      &
     &                       r2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_u2dvar(ng),                      &
     &                       u2dvar, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_v2dvar(ng),                      &
     &                       v2dvar, 2)
# if defined SSH_TIDES || defined UV_TIDES
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_rtides(ng),                      &
     &                       r2dvar+4, 3, 1, NTC(ng))
# endif
# ifdef SOLVE3D
#  ifdef VEGETATION
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_a3dvar(ng),                      &
     &                       a3dvar, 3, 1, NVEG)
#  endif
#  ifdef SEDIMENT
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_b3dvar(ng),                      &
     &                       b3dvar, 3, 1, Nbed)
#  endif
#  if defined DIAGNOSTICS_BIO && defined ECOSIM
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_l3dvar(ng),                      &
     &                       l3dvar, 3, 1, NDbands)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_l4dvar(ng),                      &
     &                       l4dvar, 4, 1, N(ng), NDbands)
#  endif
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_p3dvar(ng),                      &
     &                       p3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_r3dvar(ng),                      &
     &                       r3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_u3dvar(ng),                      &
     &                       u3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_v3dvar(ng),                      &
     &                       v3dvar, 3, 1, N(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_w3dvar(ng),                      &
     &                       w3dvar, 3, 0, N(ng))
!
        IF (Nslice.gt.0) THEN
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,  &
     &                         ioDesZ_dp_p3dvar(ng),                    &
     &                         p3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,  &
     &                         ioDesZ_dp_r3dvar(ng),                    &
     &                         r3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,  &
     &                         ioDesZ_dp_u3dvar(ng),                    &
     &                         u3dvar, 3, 1, Nslice)
          CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,  &
     &                         ioDesZ_dp_v3dvar(ng),                    &
     &                         v3dvar, 3, 1, Nslice)
        END IF
# endif
# if defined ADJUST_STFLUX && defined DISTRIBUTE
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_r2dfrc(ng),                      &
     &                       r2dvar, 3, 1, Nfrec(ng))
# endif
# ifdef ADJUST_WSTRESS
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_u2dfrc(ng),                      &
     &                       u2dvar, 3, 1, Nfrec(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_v2dfrc(ng),                      &
     &                       v2dvar, 3, 1, Nfrec(ng))
# endif
      END DO
!
!  Set I/O decomposition descriptors for special single precision
!  restart and harmonics data.
!
      DO ng=1,Ngrids
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_rubar(ng),                       &
     &                       u2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_rvbar(ng),                       &
     &                       v2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_rzeta(ng),                       &
     &                       r2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_ubar(ng),                        &
     &                       u2dvar, 3, 1, 3)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_vbar(ng),                        &
     &                       v2dvar, 3, 1, 3)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_zeta(ng),                        &
     &                       r2dvar, 3, 1, 3)
# ifdef SOLVE3D
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_ruvel(ng),                       &
     &                       u3dvar, 4, 0, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_rvvel(ng),                       &
     &                       v3dvar, 4, 0, N(ng), 1, 2)
#  if defined GLS_MIXING || defined MY25_MIXING
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_tkevar(ng),                      &
     &                       r3dvar, 4, 0, N(ng), 1, 2)
#  endif
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_trcvar(ng),                      &
     &                       r3dvar, 4, 1, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_uvel(ng),                        &
     &                       u3dvar, 4, 1, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_vvel(ng),                        &
     &                       v3dvar, 4, 1, N(ng), 1, 2)
# endif

# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_r2dhar(ng),                      &
     &                       r2dvar, 3, 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_u2dhar(ng),                      &
     &                       u2dvar, 3, 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_v2dhar(ng),                      &
     &                       v2dvar, 3, 0, 2*NTC(ng))
#  ifdef SOLVE3D
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_r3dhar(ng),                      &
     &                       r3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_u3dhar(ng),                      &
     &                       u3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_real,      &
     &                       ioDesc_sp_v3dhar(ng),                      &
     &                       v3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
#  endif
# endif
      END DO
!
!  Set I/O decomposition descriptors for special double precision
!  restart and harmonics data.
!
      DO ng=1,Ngrids
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_rubar(ng),                       &
     &                       u2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_rvbar(ng),                       &
     &                       v2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_rzeta(ng),                       &
     &                       r2dvar, 3, 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_ubar(ng),                        &
     &                       u2dvar, 3, 1, 3)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_vbar(ng),                        &
     &                       v2dvar, 3, 1, 3)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_zeta(ng),                        &
     &                       r2dvar, 3, 1, 3)
# ifdef SOLVE3D
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_ruvel(ng),                       &
     &                       u3dvar, 4, 0, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_rvvel(ng),                       &
     &                       v3dvar, 4, 0, N(ng), 1, 2)
#  if defined GLS_MIXING || defined MY25_MIXING
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_tkevar(ng),                      &
     &                       r3dvar, 4, 0, N(ng), 1, 2)
#  endif
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_trcvar(ng),                      &
     &                       r3dvar, 4, 1, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_uvel(ng),                        &
     &                       u3dvar, 4, 1, N(ng), 1, 2)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_vvel(ng),                        &
     &                       v3dvar, 4, 1, N(ng), 1, 2)
# endif

# if defined AVERAGES  && defined AVERAGES_DETIDE && \
    (defined SSH_TIDES || defined UV_TIDES)
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_r2dhar(ng),                      &
     &                       r2dvar, 3, 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_u2dhar(ng),                      &
     &                       u2dvar, 3, 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_v2dhar(ng),                      &
     &                       v2dvar, 3, 0, 2*NTC(ng))
#  ifdef SOLVE3D
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_r3dhar(ng),                      &
     &                       r3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_u3dhar(ng),                      &
     &                       u3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
        CALL field_iodecomp (ng, pioSystem(IpioROMS,ng), PIO_double,    &
     &                       ioDesc_dp_v3dhar(ng),                      &
     &                       v3dvar, 4, 1, N(ng), 0, 2*NTC(ng))
#  endif
# endif
      END DO
!
      RETURN
      END SUBROUTINE set_iodecomp

# if defined PROPAGATOR && defined CHECKPOINTING
!
      SUBROUTINE state_iodecomp (ng, ioSystem, ioType, ioDesc,          &
     &                           ioVname, ndims)
!
!***********************************************************************
!                                                                      !
!  Sets the IO decomposition descriptor for ROMS packed state variable !
!  types.                                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     ioSystem     PIO system descriptor (TYPE IOSystem_desc_t)        !
!     ioType       PIO kind variable type (integer)                    !
!     ioVname      State variable name (string)                        !
!     ndims        Number of state variable dimensions (integer)       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ioDesc       IO decomposition descriptor (TYPE io_desc_t)        !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ioType, ndims
!
      character (len=*) :: ioVname
!
      TYPE (IOSystem_desc_t), intent(in) :: ioSystem
      TYPE (io_desc_t), intent(out) :: ioDesc
!
!  Local variable declarations.
!
      integer :: Is, Ie, Isize, Js, Je, Jsize, my_size
      integer :: Imin, Imax, Ioff, Jmin, Jmax, Joff
      integer :: i, ic, j, jc, np

      integer(PIO_Offset_kind), allocatable :: map_decomp(:)
!
!-----------------------------------------------------------------------
!  Set the PIO computational decomposition for ROMS packed state
!  variables.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(ioVname))
        CASE ('Bvec')
          Is=Nstr(ng)
          Ie=Nend(ng)
          Js=1
          Je=NCV
          Ioff=0
          Joff=1
        CASE ('resid')
          Is=Nstr(ng)
          Ie=Nend(ng)
        CASE ('SworkD')
          Is=MyRank*3*Nstate(ng)+1
          Ie=MIN(Is+3*Nstate(ng)-1, 3*Mstate(ng))
      END SELECT
!
!  Starting/ending I- and J-indices for each decomposition tile.
!
      IF (ndims.eq.1) THEN
        Imin=Is
        Imax=Ie
        Isize=Ie-Is+1
        my_size=Isize
      ELSE IF (ndims.eq.2) THEN
        Imin=Is
        Imax=Ie
        Jmin=Js
        Jmax=Je
        Isize=Ie-Is+1
        Jsize=Je-Js+1
        my_size=Isize*Jsize
      END IF
!
!  Allocate 1D array for mapping of the storage order of the variable to
!  its memory order.
!
      IF (.not.ALLOCATED(map_decomp)) THEN
        allocate ( map_decomp(my_size) )
      END IF
      map_decomp=0_PIO_Offset_kind
!
!  Set variable decomposition mapping.
!
      IF (ndims.eq.1) THEN
        np=0
        DO i=Imin,Imax
          np=np+1
          map_decomp(np)=i
        END DO
      ELSE IF (ndims.eq.2) THEN
        np=0
        DO j=Jmin,Jmax
          jc=(j-Joff)*Isize
          DO i=Imin,Imax
            np=np+1
            ic=i+Ioff+jc
            map_decomp(np)=ic
          END DO
        END DO
      END IF
!
!  Set IO decomposition descriptor
!
      IF (ndims.eq.1) THEN
        CALL PIO_InitDecomp (ioSystem, ioType, (/Isize/),               &
     &                       map_decomp, ioDesc)
      ELSE IF (ndims.eq.2) THEN
        CALL PIO_InitDecomp (ioSystem, ioType, (/Isize,Jsize/),         &
     &                       map_decomp, ioDesc)
      END IF
!
!  Deallocate.
!
      IF (ALLOCATED(map_decomp)) deallocate (map_decomp)
!
      RETURN
      END SUBROUTINE state_iodecomp
# endif
#endif
      END MODULE set_pio_mod
