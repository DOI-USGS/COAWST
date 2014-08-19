#include "wrfcpp.h"
      MODULE atm_coupler_mod
#ifdef WRF_COUPLING
!svn $Id$
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between WRF    !
!  other coupled model(s) using the Model Coupling Toolkit (MCT).      !
!                                                                      !
!=======================================================================
!
!  Componenet model registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain decompositin descriptor datatype and assocoiated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field storage data types and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_indxR => indexRA
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communitcations scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY: MCT_send => send
      USE m_Transfer, ONLY: MCT_recv => recv
      USE m_Transfer, ONLY: MCT_isend => isend
      USE m_Transfer, ONLY: MCT_irecv => irecv
      USE m_Transfer, ONLY: MCT_waitr => waitrecv
      USE m_Transfer, ONLY: MCT_waits => waitsend
!
      USE mct_wrf_coupler_params
      USE module_parallel
!
      implicit none
!
      PRIVATE

      PUBLIC :: INITIALIZE_ATM_ROUTERS
      PUBLIC :: initialize_atm_coupling
      PUBLIC :: atm_coupling
# ifdef ROMS_COUPLING
      PUBLIC :: atm2ocn_coupling
      PUBLIC :: atmfocn_coupling
      PUBLIC :: atm_coupling_aux4
# endif
# ifdef SWAN_COUPLING
      PUBLIC :: atm2wav_coupling
      PUBLIC :: atmfwav_coupling
# endif
      PUBLIC :: finalize_atm_coupling

      include 'mpif.h'
!
!  Declarations.
!

      TYPE T_GlobalSegMap_G
        TYPE(GlobalSegMap) :: GSMapWRF         ! GloabalSegMap variables
      END TYPE T_GlobalSegMap_G
      TYPE (T_GlobalSegMap_G), ALLOCATABLE :: GlobalSegMap_G(:)

# if defined ROMS_COUPLING
      TYPE T_AttrVect_O
        TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
        TYPE(AttrVect) :: ocn2atm_AV 
      END TYPE T_AttrVect_O
      TYPE (T_AttrVect_O), ALLOCATABLE :: AttrVect_O(:,:)
# endif

# if defined SWAN_COUPLING
      TYPE T_AttrVect_W
       TYPE(AttrVect) :: atm2wav_AV            ! AttrVect variables
       TYPE(AttrVect) :: wav2atm_AV            ! AttrVect variables
      END TYPE T_AttrVect_W
      TYPE (T_AttrVect_W), ALLOCATABLE :: AttrVect_W(:)
# endif

# if defined SWAN_COUPLING
      TYPE T_Router_W
        type(Router)   :: WRFtoSWAN           ! Router variables
      END TYPE T_Router_W
      TYPE (T_Router_W), ALLOCATABLE :: Router_W(:,:)
# endif
# if defined ROMS_COUPLING
      TYPE T_Router_O
        type(Router)   :: WRFtoROMS           ! Router variables
      END TYPE T_Router_O
      TYPE (T_Router_O), ALLOCATABLE :: Router_O(:,:)
# endif

      CONTAINS

      SUBROUTINE initialize_atm_coupling(grid, ia)
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none
!
!  Imported variable definitions.
!
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia

      include 'mpif.h'
!
!  Local variable declarations.  
!
      integer :: MyError, MyRank
      integer :: npoints, gsmsize, nprocs, localsize
      integer :: j, jc, Isize, Jsize
      integer :: is, ie, js, je, cid, cad, io

      integer, pointer :: start(:), length(:)
      character (len=120)  :: to_add, avstring

!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (ATM_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (ATM_COMM_WORLD, nprocs, MyError)
!
      IF (ia.eq.1) THEN
        ALLOCATE(GlobalSegMap_G(Natm_grids))
# if defined ROMS_COUPLING
        ALLOCATE(AttrVect_O(Natm_grids,Nocn_grids))
# endif
# if defined SWAN_COUPLING
        ALLOCATE(AttrVect_W(Natm_grids))
# endif
      END IF
!
!  Initialize MCT coupled model registry.
!
      ATMid=atmids(ia)
      IF (Natm_grids.gt.1) THEN
        CALL MCTWorld_init (N_mctmodels,MPI_COMM_WORLD,                 &
     &                      ATM_COMM_WORLD,myids=atmids)
      ELSE
        CALL MCTWorld_init (N_mctmodels,MPI_COMM_WORLD,                 &
     &                      ATM_COMM_WORLD,ATMid)
      END IF
!
!  Initialize a Global Segment Map for non-haloed transfer of data out
!  of WRF. Determine non-haloed start and length arrays for this
!  processor.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!  Determine tile size
!
      Isize=ie-is+1
      Jsize=je-js+1
      allocate( start(Jsize) )
      allocate( length(Jsize) )
      jc=0
      DO j=js,je
        jc=jc+1
        start(jc)=(j-1)*(grid%ed31-1)+is
        length(jc)=Isize
      END DO
      gsmsize=Isize*Jsize
!
      CALL GlobalSegMap_init (GlobalSegMap_G(ia)%GSMapWRF,              &
     &                        start, length, 0, ATM_COMM_WORLD, ATMid)
# ifdef ROMS_COUPLING
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model.
      cad=LEN(avstring)
      DO j=1,cad
        avstring(j:j)=''
      END DO
      cid=1
!
      to_add='GSW'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':GLW'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
#  ifdef ATM2OCN_FLUXES
      to_add=':LH'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':HFX'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':USTRESS'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':VSTRESS'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
      to_add=':MSLP'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM || \
   (defined SHORTWAVE && defined ANA_SRFLUX)
      to_add=':RELH'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':T2'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined BULK_FLUXES || defined ECOSIM
      to_add=':U10'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
!
      to_add=':V10'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  ifdef CLOUDS
      to_add=':CLDFRA'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if !defined ANA_RAIN && defined EMINUSP
      to_add=':RAIN'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
!
#  if defined EMINUSP
      to_add=':EVAP'
      cad=LEN_TRIM(to_add)
      write(avstring(cid:cid+cad-1),'(a)') to_add(1:cad)
      cid=cid+cad
#  endif
      cad=LEN_TRIM(avstring)
      avstring=avstring(1:cad)
!
      DO io=1,Nocn_grids
        CALL AttrVect_init (AttrVect_O(ia,io)%atm2ocn_AV,               &
     &                      rlist=TRIM(avstring),lsize=gsmsize)
        CALL AttrVect_zero (AttrVect_O(ia,io)%atm2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
# ifdef MCT_INTERP_OC2AT
        CALL AttrVect_init (AttrVect_O(ia,io)%ocn2atm_AV,               &
     &                      rList="SST:CPL_MASK",lsize=gsmsize)
# else
        CALL AttrVect_init (AttrVect_O(ia,io)%ocn2atm_AV,               &
     &                      rList="SST",lsize=gsmsize)
# endif
        CALL AttrVect_zero (AttrVect_O(ia,io)%ocn2atm_AV)
      END DO
# endif
# ifdef SWAN_COUPLING
!
!  Initialize attribute vector holding the export data code string of
!  the wave model.
!
      CALL AttrVect_init (AttrVect_W(ia)%atm2wav_AV,                    &
     &                    rlist="U10:V10",lsize=gsmsize)
      CALL AttrVect_zero (AttrVect_W(ia)%atm2wav_AV)
!
!  Initialize attribute vector for data from SWAN.
!
#  ifdef MCT_INTERP_WV2AT
      CALL AttrVect_init (AttrVect_W(ia)%wav2atm_AV,                    &
     &                    rList="HSIGN:WLENP:RTP:CPL_MASK",             &
     &                    lsize=gsmsize)
#  else
      CALL AttrVect_init (AttrVect_W(ia)%wav2atm_AV,                    &
     &                    rList="HSIGN:WLENP:RTP",lsize=gsmsize)
#  endif
      CALL AttrVect_zero (AttrVect_W(ia)%wav2atm_AV)
# endif

      deallocate (start)
      deallocate (length)

      RETURN
      END SUBROUTINE initialize_atm_coupling

      SUBROUTINE INITIALIZE_ATM_ROUTERS
!
!=======================================================================
!                                                                      !
!  Initialize atm routers.                                             !
!                                                                      !
!=======================================================================

      implicit none
!
!  Local variable declarations.
!
      integer :: io, iw, ia
!
!  Initialize MCT Routers.
!
# ifdef ROMS_COUPLING
      ALLOCATE(Router_O(Natm_grids,Nocn_grids))
!
!  Initialize a router to the ocean model component.
!
      DO io=1,Nocn_grids
        DO ia=1,Natm_grids
          OCNid=ocnids(io)
          CALL Router_init (OCNid, GlobalSegMap_G(ia)%GSMapWRF,         &
     &                      ATM_COMM_WORLD, Router_O(ia,io)%WRFtoROMS)
        END DO
      END DO
# endif
# ifdef SWAN_COUPLING
      ALLOCATE(Router_W(Natm_grids,Nwav_grids))
!
!  Initialize a router to the wav model component.
!
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          WAVid=wavids(iw)
          CALL Router_init (WAVid, GlobalSegMap_G(ia)%GSMapWRF,         &
     &                      ATM_COMM_WORLD, Router_W(ia,iw)%WRFtoSWAN)
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE INITIALIZE_ATM_ROUTERS

      SUBROUTINE atm_coupling(grid, num_steps)
!
!=======================================================================
!                                                                      !
!  Determine which other model to do the coupling at this time step.   !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none

!
!  Imported variable definitions.
!
      integer, intent(in) :: num_steps
!     TYPE(domain) , INTENT (IN) :: grid 
      TYPE(domain) , POINTER     :: grid 
      TYPE(domain) , POINTER     :: grid_ptr
!
!  Local variable declarations.  
!
      integer :: io, iw, ia, offset
      integer :: kid, num_ksteps
!
      IF (num_steps.eq.0) THEN
        offset=0
      ELSE
!       offset=-1
        offset=0
      END IF
# ifdef ROMS_COUPLING
      IF (MOD(num_steps+offset, nATM2OCN(1,1)).eq.0) THEN
        DO io=1,Nocn_grids
          ia=1
          CALL atm2ocn_coupling (grid,ia,io)
          grid_ptr => grid
          DO WHILE ( ASSOCIATED( grid_ptr ) )
            DO kid = 1, max_nests
              IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                ia=ia+1
                num_ksteps=grid_ptr%itimestep
                CALL atm2ocn_coupling (grid_ptr%nests(kid)%ptr,ia,io)
              END IF
            END DO
            grid_ptr => grid_ptr%sibling
          END DO
        END DO
      END IF
!
      IF (MOD(num_steps+offset, nATMFOCN(1,1)).eq.0) THEN
        DO io=1,Nocn_grids
          ia=1
          CALL atmfocn_coupling (grid,ia,io,0)
          grid_ptr => grid
          DO WHILE ( ASSOCIATED( grid_ptr ) )
            DO kid = 1, max_nests
              IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                ia=ia+1
                num_ksteps=grid_ptr%itimestep
                CALL atmfocn_coupling (grid_ptr%nests(kid)%ptr,ia,io,0)
              END IF
            END DO
            grid_ptr => grid_ptr%sibling
          END DO
        END DO
      END IF
# endif
# ifdef SWAN_COUPLING
      IF (MOD(num_steps+offset, nATM2WAV(1,1)).eq.0) THEN
        DO iw=1,Nwav_grids
          ia=1
          CALL atm2wav_coupling (grid,ia,iw)
          grid_ptr => grid
          DO WHILE ( ASSOCIATED( grid_ptr ) )
            DO kid = 1, max_nests
              IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                ia=ia+1
                num_ksteps=grid_ptr%itimestep
                CALL atm2wav_coupling (grid_ptr%nests(kid)%ptr,ia,iw)
              END IF
            END DO
            grid_ptr => grid_ptr%sibling
          END DO
        END DO
      END IF
!
      IF (MOD(num_steps+offset, nATMFWAV(1,1)).eq.0) THEN
        DO iw=1,Nwav_grids
          ia=1
          CALL atmfwav_coupling (grid,ia,iw)
          grid_ptr => grid
          DO WHILE ( ASSOCIATED( grid_ptr ) )
            DO kid = 1, max_nests
              IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
                ia=ia+1
                num_ksteps=grid_ptr%itimestep
                CALL atmfwav_coupling (grid_ptr%nests(kid)%ptr,ia,iw)
              END IF
            END DO
            grid_ptr => grid_ptr%sibling
          END DO
        END DO
      END IF
# endif

      RETURN
      END SUBROUTINE atm_coupling

      SUBROUTINE atm_coupling_aux4(grid, num_steps)
!
!=======================================================================
!                                                                      !
!  Determine which other model to do the coupling at this time step.   !
!                                                                      !
!=======================================================================
!
      USE module_domain

      implicit none
!
!  Imported variable definitions.
!
      integer, intent(in) :: num_steps
      TYPE(domain) , POINTER     :: grid 
      TYPE(domain) , POINTER     :: grid_ptr
!
!  Local variable declarations.  
!
      integer :: io, ia
      integer :: kid, num_ksteps
!
# ifdef ROMS_COUPLING
      DO io=1,Nocn_grids
        ia=1
        CALL atmfocn_coupling (grid,ia,io,1)
        grid_ptr => grid
        DO WHILE ( ASSOCIATED( grid_ptr ) )
          DO kid = 1, max_nests
            IF ( ASSOCIATED( grid_ptr%nests(kid)%ptr ) ) THEN
              ia=ia+1
              CALL atmfocn_coupling (grid_ptr%nests(kid)%ptr,ia,io,1)
            END IF
          END DO
          grid_ptr => grid_ptr%sibling
        END DO
      END DO
# endif

      RETURN
      END SUBROUTINE atm_coupling_aux4

# ifdef ROMS_COUPLING
      SUBROUTINE atm2ocn_coupling (grid, ia, io)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, io
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET

      real, parameter :: eps=1.0e-10
      real, pointer :: AA(:)
      real :: cff1, cff2, cff3, rnum, rden, c04, c05
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
!
!-----------------------------------------------------------------------
!  Export fields from atmosphere (WRF) to ocean (ROMS) model.
!-----------------------------------------------------------------------
!     GSW        Short wave raditaion (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%GSW(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "GSW",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     GLW        Long wave raditaion (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%GLW(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "GLW",   &
     &                           AA, Asize)
#  ifdef ATM2OCN_FLUXES
!-----------------------------------------------------------------------
!     LH     Latent heat flux (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%LH(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "LH",    &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     HFX     Sensible heat flux (W m-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%HFX(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "HFX",   &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     USTRESS    Surface u-stress (m2 s-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          cff1=1.0/(grid%alt(i,1,j)+eps)
          cff2=2.0/(((grid%u_2(i,1,j)+grid%u_2(i+1,1,j))**2+            &
     &               (grid%v_2(i,1,j)+grid%v_2(i,1,j+1))**2)**0.5+eps)
          ij=ij+1
          AA(ij)=cff1*cff2*(grid%UST(i,j)**2)*                          &
     &           0.5*(grid%u_2(i,1,j)+grid%u_2(i+1,1,j))
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV,"USTRESS",&
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     VSTRESS    Surface v-stress (m2 s-2).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          cff1=1.0/(grid%alt(i,1,j)+eps)
          cff2=2.0/(((grid%u_2(i,1,j)+grid%u_2(i+1,1,j))**2+            &
     &               (grid%v_2(i,1,j)+grid%v_2(i,1,j+1))**2)**0.5+eps)
          ij=ij+1
          AA(ij)=cff1*cff2*(grid%UST(i,j)**2)*                          &
     &           0.5*(grid%v_2(i,1,j)+grid%v_2(i,1,j+1))
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV,"VSTRESS",&
     &                           AA, Asize)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!-----------------------------------------------------------------------
!     MSLP       Surface atmospheric pressure (Pa).
!     Use the hypsometric equation to reduce 
!     surface pressure to mean sea level pressure.
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%PSFC(i,j)*                                        &
     &           exp((9.81*grid%ht(i,j))/                               &
     &           (287.0*grid%T2(i,j)*(1.0+0.61*grid%Q2(i,j))))
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "MSLP",  &
     &                           AA, Asize)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || \
     (defined SHORTWAVE && defined ANA_SRFLUX)
!-----------------------------------------------------------------------
!     RELH       Surface air relative humidity (-).
!
      ij=0
      DO j=js,je
        DO i=is,ie
!
!         Calculate 2-m pressure using hypsometric equation. 
!         Assume temp at 2m = temp at 0m.
!
          cff1 = grid%PSFC(i,j) / (exp((9.81*2.0)/(287.0*grid%T2(i,j))))
!
!         Compute specific humidity using the 2-m mixing ratio and 2-m pressure.
!
          rnum = grid%Q2(i,j)*cff1
          rden  = (grid%Q2(i,j)*(1.-0.622)+0.622)
          cff2 = rnum/rden                         ! e is the specific humidity
!
!         Compute saturation specific humidity using Bolton equation 10.
!
          c04 = 17.67*(grid%T2(i,j)-273.15)
          c05 = (grid%T2(i,j)-273.15) + 243.5
          cff3  = 6.112*exp(c04/c05)
!
          ij=ij+1
          AA(ij)=cff2/cff3
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "RELH",  &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     T2         Surface 2m air temperature (Convert to C).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%T2(i,j)-273.15
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "T2",    &
     &                           AA, Asize)
#  endif
!-----------------------------------------------------------------------
#  if defined BULK_FLUXES || defined ECOSIM
!     U10        U-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%U10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "U10",   &
     &                           AA, Asize)
#  endif
!-----------------------------------------------------------------------
#  if defined BULK_FLUXES || defined ECOSIM
!     V10        V-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%V10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "V10",   &
     &                           AA, Asize) 
#  endif
#  ifdef CLOUDS
!-----------------------------------------------------------------------
!     CLDFRA     Cloud fraction (--, 0-1.0).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%CLDFRA(i,1,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "CLDFRA",&
     &                           AA, Asize)
#  endif
#  if !defined ANA_RAIN && defined EMINUSP
!-----------------------------------------------------------------------
!     RAIN       Precipitation (Convert to m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=0.001*(grid%RAINCV(i,j)+grid%RAINNCV(i,j))/grid%dt
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "RAIN",  &
     &                           AA, Asize)
#  endif
#  if defined EMINUSP
!-----------------------------------------------------------------------
!     EVAP      Evaporation (kg/m2 to Convert to m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
!         AA(ij)=0.001*(grid%SFCEVP(i,j))/grid%dt
          AA(ij)=0.001*(grid%QFX(i,j))
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_O(ia,io)%atm2ocn_AV, "EVAP",  &
     &                           AA, Asize)
#  endif
!-----------------------------------------------------------------------
!  Send fields to ocean model.
!
      Tag=io*100+ia*10+0
      CALL MCT_ISEND (AttrVect_O(ia,io)%atm2ocn_AV,                     &
     &                Router_O(ia,io)%WRFtoROMS, Tag)
      CALL MCT_waits (Router_O(ia,io)%WRFtoROMS)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,36) ' ## WRF grid ',ia,                                &
     &                    ' sent data to ROMS grid ',io
 36     FORMAT (a14,i2,a24,i2)
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_atm_coupling
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      RETURN
      END SUBROUTINE atm2ocn_coupling
# endif

# ifdef SWAN_COUPLING
      SUBROUTINE atm2wav_coupling (grid, ia, iw)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, iw
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET

      real, parameter :: eps=1.0e-10
      real, pointer :: AA(:)
      real :: cff1, cff2, cff3, rnum, rden, c04, c05
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
!-----------------------------------------------------------------------
!  Send U10 and V10 to SWAN.
!-----------------------------------------------------------------------
!     U10        U-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%U10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_W(ia)%atm2wav_AV, "U10",      &
     &                           AA, Asize)
!-----------------------------------------------------------------------
!     V10        V-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          AA(ij)=grid%V10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (AttrVect_W(ia)%atm2wav_AV, "V10",      &
     &                           AA, Asize)
!
      Tag=0+ia*10+iw
      CALL MCT_ISEND (AttrVect_W(ia)%atm2wav_AV,                        &
     &                Router_W(ia,iw)%WRFtoSWAN, Tag)
      CALL MCT_waits (Router_W(ia,iw)%WRFtoSWAN)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,37) ' ## WRF grid ',ia,                                &
     &                    ' sent data to SWAN grid ',iw
 37     FORMAT (a14,i2,a24,i2)
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_atm_coupling
      END IF
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      RETURN
      END SUBROUTINE atm2wav_coupling
# endif

# ifdef ROMS_COUPLING
      SUBROUTINE atmfocn_coupling (grid, ia, io, aux4flag)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, io, aux4flag
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET

      real, parameter :: eps=1.0e-10
      real, pointer :: AA(:)
      real, pointer :: Amask(:)
      real :: cff1, cff2, cff3, rnum, rden, c04, c05
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
      allocate ( Amask(Asize), stat=ierr )
!
!  Schedule receiving fields from ocean model.
!  For now, use this aux4flag to do update if sst_update is used.
!
      IF (aux4flag.eq.0) THEN
        Tag=io*100+ia*10+0
        CALL MCT_irecv (AttrVect_O(ia,io)%ocn2atm_AV,                   &
     &                  Router_O(ia,io)%WRFtoROMS,Tag)
!
!       Wait to make sure the OCN data has arrived.
!
        CALL MCT_waitr (AttrVect_O(ia,io)%ocn2atm_AV,                   &
     &                  Router_O(ia,io)%WRFtoROMS)
!
        IF (MYRANK.EQ.0) THEN
          WRITE (*,38) ' ## WRF grid ',ia,                              &
     &                      ' recvd data from ROMS grid ',io
 38       FORMAT (a13,i2,a27,i2)
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling fail wrfcplr, MyStatus= ', MyError
          CALL finalize_atm_coupling
        END IF
      END IF
!
!  SST (Convert to K).
!  USE CPLMASK:  coupling mask (0 for data read in wrflowinput,
!                               1 data received from the coupler)
#  ifdef MCT_INTERP_OC2AT
      CALL AttrVect_exportRAttr (AttrVect_O(ia,io)%ocn2atm_AV,          &
     &                           "CPL_MASK", AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#    ifdef SST_CONST
          Amask(ij)=0.
#    else
          Amask(ij)=AA(ij)
#    endif
        END DO
      END DO
#  else
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#    ifdef SST_CONST
          Amask(ij)=0.
#    else
          Amask(ij)=1.
#    endif
        END DO
      END DO
#  endif
!
      CALL AttrVect_exportRAttr (AttrVect_O(ia,io)%ocn2atm_AV,          &
     &                           "SST", AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          grid%sst(i,j)=(AA(ij)+273.15)*Amask(ij)+                      &
     &                    grid%sst(i,j)*(1.0-Amask(ij))
        END DO
      END DO
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      deallocate (Amask)
      RETURN
      END SUBROUTINE atmfocn_coupling
# endif

# ifdef SWAN_COUPLING
      SUBROUTINE atmfwav_coupling (grid, ia, iw)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupled data streams.          !
!  Currently, the following data streams are processed:                !
!                                                                      !
!  Possible fields exported to the OCEAN Model:                        !
!                                                                      !
!     * GSW        Short wave raditaion  (Watts/m2)                    !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * MSLP       Mean Sea Level Pressure (Pa)                        !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * EVAP       Evaporation          (m/s)                          !
!                                                                      !
!  Fields exported to the WAVE Model:                                  !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!                                                                      !
!  Fields acquired from the WAVE Model:                                !
!                                                                      !
!     * HISGN      Significant wave heigth (m)                         !
!     * WLENP      Peak wave length (m)                                !
!     * RTP        Peak wave period (s)                                !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
      integer, intent(in) :: ia, iw
!
!  Local variable declarations.
!
      integer :: is, ie, js, je, ij, Tag
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET

      real, parameter :: eps=1.0e-10
      real, pointer :: AA(:)
#  ifdef MCT_INTERP_WV2AT
      real, pointer :: Amask(:)
#  endif
      real :: cff1, cff2, cff3, rnum, rden, c04, c05
!  Set grid range.
!
      is = grid%sp31
      ie = grid%ep31
      js = grid%sp33
      je = grid%ep33
      IF (grid%ed31.eq.ie) THEN
        ie=ie-1
      END IF
      IF (grid%ed33.eq.je) THEN
        je=je-1
      END IF
!
!-----------------------------------------------------------------------
!  Send atmosphere fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (ATM_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (ATM_COMM_WORLD, nprocs, MyError)
!
!  Get the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GlobalSegMap_G(ia)%GSMapWRF,             &
     &                         ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
#  ifdef MCT_INTERP_WV2AT
      allocate ( Amask(Asize), stat=ierr )
#  endif
!
!  Schedule receiving fields from wave model.
!
      Tag=0+ia*10+iw
      CALL MCT_irecv (AttrVect_W(ia)%wav2atm_AV,                        &
     &                Router_W(ia,iw)%WRFtoSWAN, Tag)
!
!     Wait to make sure the WAV data has arrived.
!
      CALL MCT_waitr (AttrVect_W(ia)%wav2atm_AV,                        &
     &                 Router_W(ia,iw)%WRFtoSWAN)
!
      IF (MYRANK.EQ.0) THEN
        WRITE (*,40) ' ## WRF grid ',ia,                                &
     &                    ' recvd  data from SWAN grid ',iw
 40     FORMAT (a14,i2,a24,i2)
      END IF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupl fail wrfcplr, MyStatus= ', MyError
        CALL finalize_atm_coupling
      END IF
!
!  CPL_MASK
!
#  ifdef MCT_INTERP_WV2AT
      CALL AttrVect_exportRAttr (AttrVect_W(ia)%wav2atm_AV, "CPL_MASK", &
     &                           AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
          Amask(ij)=AA(ij)
        END DO
      END DO
#  endif
!
!  Hwave.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia)%wav2atm_AV, "HSIGN",    &
     &                           AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#  ifdef MCT_INTERP_WV2AT
          grid%hwave(i,j)=MAX(0.0,AA(ij))*Amask(ij)+                    &
     &                    grid%hwave(i,j)*(1.0-Amask(ij))
#  else
          grid%hwave(i,j)=MAX(0.0,AA(ij))
#  endif
        END DO
      END DO
!
!  Lwave.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia)%wav2atm_AV, "WLENP",    &
     &                           AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#  ifdef MCT_INTERP_WV2AT
          grid%lwavep(i,j)=MAX(0.0,AA(ij))*Amask(ij)+                   &
     &                     grid%lwavep(i,j)*(1.0-Amask(ij))
#  else
          grid%lwavep(i,j)=MAX(0.0,AA(ij))
#  endif
        END DO
      END DO
!
!  Pwave_top.
!
      CALL AttrVect_exportRAttr (AttrVect_W(ia)%wav2atm_AV, "RTP",      &
     &                           AA, Asize)
      ij=0
      DO j=js,je
        DO i=is,ie
          ij=ij+1
#  ifdef MCT_INTERP_WV2AT
          grid%pwave(i,j)=MAX(1.0,AA(ij))*Amask(ij)+                    &
     &                     grid%pwave(i,j)*(1.0-Amask(ij))
#  else
          grid%pwave(i,j)=MAX(1.0,AA(ij))
#  endif
        END DO
      END DO
!
!  Deallocate communication arrays.
!
      deallocate (AA)
#  ifdef MCT_INTERP_WV2AT
      deallocate (Amask)
#  endif
      RETURN
      END SUBROUTINE atmfwav_coupling
#endif

      SUBROUTINE finalize_atm_coupling
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ia, io, MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
!     CALL Router_clean (RoutWRFtoROMS)
      DO ia=1,Natm_grids
        CALL GlobalSegMap_clean (GlobalSegMap_G(ia)%GSMapWRF)
# ifdef ROMS_COUPLING
        DO io=1,Nocn_grids
          CALL AttrVect_clean (AttrVect_O(ia,io)%atm2ocn_AV)
          CALL AttrVect_clean (AttrVect_O(ia,io)%ocn2atm_AV)
        END DO
# endif
# ifdef SWAN_COUPLING
        CALL AttrVect_clean (AttrVect_W(ia)%atm2wav_AV)
        CALL AttrVect_clean (AttrVect_W(ia)%wav2atm_AV)
# endif
      END DO
!     CALL MCTWorld_clean ()

      END SUBROUTINE finalize_atm_coupling
#endif
      END MODULE atm_coupler_mod
