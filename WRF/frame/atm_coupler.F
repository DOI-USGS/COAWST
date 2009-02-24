#include "cppdefs.wrf"
      MODULE atm_coupler_mod
#ifdef ATM_OCEAN
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
      USE m_Transfer, ONLY : MCT_Send => send
      USE m_Transfer, ONLY : MCT_Recv => recv
!
      implicit none
!
      PRIVATE

      PUBLIC :: initialize_ocean_coupling
      PUBLIC :: atm2ocn_coupling
      PUBLIC :: finalize_ocean_coupling

      include 'mpif.h'
!
!  Declarations.
!
      TYPE(GlobalSegMap) :: GSMapWRF          ! GloabalSegMap variables
 
      TYPE(AttrVect) :: atm2ocn_AV             ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV 
      type(Router) :: RoutWRFtoROMS           ! Router variables
# if defined SWAN_COUPLING
      TYPE(AttrVect) :: atm2wav_AV             ! AttrVect variables
      type(Router) :: RoutWRFtoSWAN           ! Router variables
# endif


# if defined SWAN_COUPLING
      integer :: ncomps = 3
# else
      integer :: ncomps = 2
# endif
      integer :: OCNid = 1
      integer :: ATMid = 2
      integer :: WAVid = 3

      CONTAINS

      SUBROUTINE initialize_ocean_coupling(grid)
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
!     USE module_dm
      USE module_parallel

      implicit none

      TYPE(domain) , INTENT (IN) :: grid 

      include 'mpif.h'
!
!  Local variable declarations.  
!
      integer :: MyError, MyRank
      integer :: npoints, gsmsize, nprocs, localsize
      integer :: j, jc, Isize, Jsize
      integer :: ids, ide, jds, jde

      integer, pointer :: start(:), length(:)
      character (len=200)   :: avstring

!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
      CALL mpi_comm_rank (ATM_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (ATM_COMM_WORLD, nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
      CALL MCTWorld_init (ncomps,MPI_COMM_WORLD,ATM_COMM_WORLD,AtmId)
!
!  Initialize a Global Segment Map for non-haloed transfer of data out
!  of WRF. Determine non-haloed start and length arrays for this
!  processor.
!
      ids = grid%sp31
      ide = grid%ep31
      jds = grid%sp33
      jde = grid%ep33
!
!  Determine tile size
!
      Isize=ide-ids+1
      Jsize=jde-jds+1
      allocate( start(Jsize) )
      allocate( length(Jsize) )
      jc=0
      DO j=jds,jde
        jc=jc+1
        start(jc)=(j-1)*grid%ed31+ids
        length(jc)=Isize
      END DO
      gsmsize=Isize*Jsize
!
      CALL GlobalSegMap_init (GSMapWRF, start, length, 0,               &
     &                        ATM_COMM_WORLD, ATMid)
!
!  Initialize attribute vector holding the export data code strings of
!  the atmosphere model.
      avstring='PSFC'
      avstring=avstring//':RELH'
      avstring=avstring//':T2'
      avstring=avstring//':U10'
      avstring=avstring//':V10'
      avstring=avstring//':CLDFRA'
      avstring=avstring//':RAIN'
      avstring=avstring//':SWDOWN'
      avstring=avstring//':GLW'
      avstring=avstring//':USTRESS'
      avstring=avstring//':VSTRESS'
      avstring=avstring//':LH'
      avstring=avstring//':HFX'
      CALL AttrVect_init (atm2ocn_AV,rlist=TRIM(avstring),lsize=gsmsize)
      CALL AttrVect_zero (atm2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
      CALL AttrVect_init (ocn2atm_AV,rList="SST",lsize=gsmsize)
      CALL AttrVect_zero (ocn2atm_AV)
!
!  Initialize a router to the ocean model component.
!
      CALL Router_init (OcnId,GSMapWRF,ATM_COMM_WORLD,RoutWRFtoROMS)

# ifdef SWAN_COUPLING
!
!  Initialize attribute vector holding the export data code string of
!  the wave model.
!
      CALL AttrVect_init (atm2wav_AV, rlist="UWind:VWind",lsize=gsmsize)
      CALL AttrVect_zero (atm2wav_AV)
!
! Initialize a router to the SWAN component.
!
      CALL Router_init (WavId,GSMapWRF,ATM_COMM_WORLD,RoutWRFtoSWAN)
# endif

      deallocate (start)
      deallocate (length)

      RETURN
      END SUBROUTINE initialize_atm2ocn_coupling

      SUBROUTINE atm2ocn_coupling (grid)
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Sea surface potential temperature (Celsius)                    !
!                                                                      !
!     * PSFC       Surface atmospheric pressure (Pa)                   !
!     * RELH       Surface air relative humidity (percent)             !
!     * T2         Surface 2m air temperature (Celsius)                !
!     * U10        U-Wind speed at 10 m (m/s)                          !
!     * V10        V-Wind speed at 10 m (m/s)                          !
!     * CLDFRA     Cloud fraction       (percent/100)                  !
!     * RAIN       Precipitation        (m/s)                          !
!     * SWDOWN     Short wave radiation (Watts/m2)                     !
!     * GLW        Long wave raditaion  (Watts/m2)                     !
!     * USTRESS    Surface U-wind stress (Pa)                          !
!     * VSTRESS    Surface v-stress      (Pa)                          !
!     * LH         Latent heat flux     (Watts/m2)                     !
!     * HFX        Sensible heat flux   (Watts/m2)                     !
!     *            Net surface heat flux (Watts/m2)                    !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * SST        Sea surface temperature                             !
!=======================================================================
!
      USE module_domain
      USE module_parallel
!
      implicit none
      TYPE(domain) , INTENT (IN) :: grid 
!
!  Local variable declarations.
!
      integer :: ids, ide, jds, jde, ij
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      integer :: Isize, Jsize, INDXG, NPROCS, OFFSET

      integer, dimension(MPI_STATUS_SIZE,4) :: status

      real, parameter :: eps=1.0e-10
      real(r8), pointer :: AA(:)
      real :: cff1, cff2, cff3, rnum, rden, c04, c05
!  Set grid range.
!
      ids = grid%sp31
      ide = grid%ep31
      jds = grid%sp33
      jde = grid%ep33
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
      Asize=GlobalSegMap_lsize(GSMapWRF,ATM_COMM_WORLD)
!
!  Allocate attribute vector array used to export/import data.
!
      allocate ( AA(Asize), stat=ierr )
!
!-----------------------------------------------------------------------
!  Export fields from atmosphere (WRF) to ocean (ROMS) model.
!-----------------------------------------------------------------------
!     PSFC       Surface atmospheric pressure (Pa).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%PSFC(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "PSFC", AA, Asize)
!
!-----------------------------------------------------------------------
!     RELH       Surface air relative humidity (-).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
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
      CALL AttrVect_importRAttr (atm2ocn_AV, "RELH", AA, Asize)
!-----------------------------------------------------------------------
!     T2         Surface 2m air temperature (Convert to C).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%T2(i,j)-273.15
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "T2", AA, Asize)
!-----------------------------------------------------------------------
!     U10        U-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%U10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "U10", AA, Asize)
# if defined SWAN_COUPLING
      CALL AttrVect_importRAttr (atm2wav_AV, "U10", AA, Asize)
# endif
!-----------------------------------------------------------------------
!     V10        V-Wind speed at 10 m (m s-1).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%V10(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "V10", AA, Asize) 
# if defined SWAN_COUPLING
      CALL AttrVect_importRAttr (atm2wav_AV, "V10", AA, Asize)
# endif
!-----------------------------------------------------------------------
!     CLDFRA     Cloud fraction (--, 0-1.0).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%CLDFRA(i,1,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "CLDFRA", AA, Asize)
!-----------------------------------------------------------------------
!     RAIN       Precipitation (Convert to m s-1).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=0.001*(grid%RAINCV(i,j)+grid%RAINNCV(i,j))
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "RAIN", AA, Asize)
!-----------------------------------------------------------------------
!     SWDOWN     Short wave radiation (W m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%SWDOWN(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "SWDOWN", AA, Asize)
!-----------------------------------------------------------------------
!     GLW        Long wave raditaion (W m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%GLW(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "GLW", AA, Asize)
!-----------------------------------------------------------------------
!     USTRESS    Surface u-stress (N m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          cff1=1.0/(grid%em_alt(i,1,j)+eps)
          cff2=1.0/((grid%em_u_2(i,1,j)**2+grid%em_v_2(i,1,j)**2)**0.5+ &
     &               eps)
          ij=ij+1
          AA(ij)=cff1*cff2*(grid%UST(i,j)**2)*(grid%em_u_2(i,1,j))
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "USTRESS", AA, Asize)
!-----------------------------------------------------------------------
!     VSTRESS    Surface v-stress (N m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          cff1=1.0/(grid%em_alt(i,1,j)+eps)
          cff2=1.0/((grid%em_u_2(i,1,j)**2+grid%em_v_2(i,1,j)**2)**0.5+ &
     &               eps)
          ij=ij+1
          AA(ij)=cff1*cff2*(grid%UST(i,j)**2)*(grid%em_v_2(i,1,j))
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "VSTRESS", AA, Asize)
!-----------------------------------------------------------------------
!     LH     Latent heat flux (W m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%LH(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "LH", AA, Asize)
!-----------------------------------------------------------------------
!     HFX     Sensible heat flux (W m-2).
!
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          AA(ij)=grid%HFX(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (atm2ocn_AV, "HFX", AA, Asize)
!-----------------------------------------------------------------------
!  Send fields to ocean model.
!
      CALL MCT_SEND (atm2ocn_AV,RoutWRFtoROMS,MyError)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,*) '== WRF sent atm fields to ROMS and Myerror= ',     &
     &              MyError
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_ocean_coupling ("Coupling failed atm_coupler")
      END IF
!
!  Schedule receiving fields from ocean model.
!
        CALL MCT_Recv (ocn2atm_AV,RoutWRFtoROMS,MyError)
        IF (MYRANK.EQ.0) THEN
          WRITE (*,*) '== WRF recvd ocean fields and Myerror= ',        &
     &                MyError
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling fail wrfcplr, MyStatus= ', MyError
          CALL finalize_ocean_coupling ("Coupling failed wrfcplr")
        END IF
!
!  SST (Convert to K).
!
      CALL AttrVect_exportRAttr (ocn2atm_AV, "SST", AA, Asize)
      ij=0
      DO j=jds,jde
        DO i=ids,ide
          ij=ij+1
          IF (AA(ij).gt.0.10) THEN
              grid%sst(i,j)=(AA(ij)+273.15)*(1-grid%landmask(i,j))+     &
     &                       grid%sst(i,j)*grid%landmask(i,j)
          END IF
        END DO
      END DO
!-----------------------------------------------------------------------
!  Send atm fields to SWAN.
!-----------------------------------------------------------------------
!
# if defined SWAN_COUPLING
      CALL MCT_SEND (atm2wav_AV,RoutWRFtoSWAN,MyError)
      IF (MYRANK.EQ.0) THEN
        WRITE (*,*) '== WRF sent atm fields to SWAN and Myerror= ',     &
     &               MyError
      ENDIF
      IF (MyError.ne.0) THEN
        WRITE (*,*) 'coupling send fail atm_coupler, error= ', MyError
        CALL finalize_ocean_coupling ("Coupling failed atm_coupler")
      END IF
# endif
!
!  Deallocate communication arrays.
!
      deallocate (AA)
      RETURN
      END SUBROUTINE atm2ocn_coupling

      SUBROUTINE finalize_atm2ocn_coupling
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
      integer :: MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (RoutWRFtoROMS)
      CALL AttrVect_clean (atm2ocn_AV)
      CALL AttrVect_clean (ocn2atm_AV)
      CALL GlobalSegMap_clean (GSMapWRF)
      CALL MCTWorld_clean ()

      END SUBROUTINE finalize_atm2ocn_coupling
#endif
      END MODULE atm_coupler_mod
