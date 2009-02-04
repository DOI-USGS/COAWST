#include "cppdefs.wrf"
      MODULE atm_coupler_mod

#if defined MODEL_COUPLING && defined MCT_LIB
!
!svn $Id: atm_coupler.F 601 2008-04-06 22:17:22Z arango $
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

      PUBLIC :: initialize_atm2ocn_coupling
      PUBLIC :: atm2ocn_coupling
      PUBLIC :: finalize_atm2ocn_coupling

      include 'mpif.h'
!
!  Declarations.
!
      TYPE(GlobalSegMap) :: GSMapWRF          ! GloabalSegMap variables
 
      TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV 
      type(Router)   :: WRFtoROMS             ! Router variables

# if defined SWAN_COUPLING
      TYPE(AttrVect) :: atm2wav_AV            ! AttrVect variables
      type(Router)   :: WRFtoSWAN             ! Router variables
# endif

      CONTAINS

      SUBROUTINE initialize_atm2ocn_coupling (grid)
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
      USE mod_coupler
!
!  Imported variable definitions.
!
      TYPE(domain), INTENT (IN) :: grid 
!
      include 'mpif.h'
!
!  Local variable declarations.
!
      integer :: MyError, MyRank
      integer :: npoints, gsmsize, nprocs, localsize
      integer :: j, jc, Isize, Jsize
      integer :: ids, ide, jds, jde

      integer, pointer :: start(:), length(:)
!
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
      CALL MCTWorld_init (Nmodels, MPI_COMM_WORLD, ATM_COMM_WORLD,      &
     &                    ATMid)
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
!       start(jc)=Isize*(j-1)+1
!       start(jc)=j*grid%ed31+ids+1
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
!
      CALL AttrVect_init (atm2ocn_AV, rlist=TRIM(ExportList(Iatmos)),   &
     &                    lsize=gsmsize)
      CALL AttrVect_zero (atm2ocn_AV)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      CALL AttrVect_init (ocn2atm_AV, rList=TRIM(ExportList(Iocean)),   &
     &                    lsize=gsmsize)
      CALL AttrVect_zero (ocn2atm_AV)
!
!  Initialize a router to the ocean model component.
!
      CALL Router_init (OCNid, GSMapWRF, ATM_COMM_WORLD, WRFtoROMS)

# ifdef SWAN_COUPLING
!
!  Initialize attribute vector holding the export data code string of
!  the wave model.
!
      CALL AttrVect_init (atm2wav_AV, rList=TRIM(ExportList(Iwaves)),   &
     &                    lsize=gsmsize)
      CALL AttrVect_zero (atm2wav_AV)
!
! Initialize a router to the SWAN component.
!
      CALL Router_init (WAVid, GSMapWRF, ATM_COMM_WORLD, WRFtoSWAN)
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
!  Fields exported to the OCEAN Model:                                 !
!                                                                      !
!     * Surface atmospheric pressure (Pa)                              !
!     * Surface air relative humidity (percent)                        !
!     * Surface (2 m) air temperature (Celsius)                        !
!     * Surface (10 m) U-wind speed (m/s)                              !
!     * Surface (10 m) V-wind speed (m/s)                              !
!     * Cloud fraction (percent/100)                                   !
!     * Precipitation (m/s)                                            !
!     * Shortwave radiation (Watts/m2)                                 !
!     * Long wave raditaion (Watts/m2)                                 !
!     * Latent heat flux (Watts/m2)                                    !
!     * Sensible heat flux (Watts/m2)                                  !
!     * Net surface heat flux (Watts/2)                                !
!     * Surface U-wind stress (Pa)                                     !
!     * Surface V-wind stress (Pa)                                     !
!                                                                      !
!  Fields imported from the OCEAN Model:                               !
!                                                                      !
!     * Sea surface potential temperature (Celsius)                    !
!                                                                      !
!=======================================================================
!
      USE module_domain
!     USE module_dm
!     USE module_io_domain
!     USE module_configure
      USE module_parallel
      USE mod_coupler
!
      implicit none
!
!  Imported variable declarations.
!
      type(domain) , intent(in) :: grid 
!
!  Local variable declarations.
!
      integer :: Asize, i, j, ids, ide, jds, jde, ij
      integer :: nprocs, MyError, MyRank, ierr
      integer :: Iexport, Iimport, id, ifield

!     integer, pointer :: points(:)

      real, parameter :: eps=1.0e-10

      real :: cff1, cff2, cff3, rnum, rden, c04, c05

      real(r8) :: RecvTime, SendTime, wtime(2)
      real(r8) :: inpbuffer(2), outbuffer(2)

      real(r8), pointer :: AA(:)

      character (len=3 ), dimension(2) :: op_handle
      character (len=20) :: wrf_time
      character (len=40) :: code
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
!     allocate ( points(Asize), stat=ierr )
!
!  Ask for points in this tile.
!
!      CALL GlobalSegMap_Ordpnts (GSMapWRF,MyRank,points)
!
!-----------------------------------------------------------------------
!  Export fields from atmosphere (WRF) to ocean (ROMS) model.
!-----------------------------------------------------------------------
!
!  Set grid range.
!
!!    ids = grid%sm31+6
!!    ide = grid%em31-6
!!    jds = grid%sm33+6
!!    jde = grid%em33-6
      ids = grid%sp31
      ide = grid%ep31
      jds = grid%sp33
      jde = grid%ep33
!
!  Initialize coupling wait time clocks.
!
      RecvTime=0.0_r8
      SendTime=0.0_r8
!
!  Schedule sending fields to the ocean model.
!
      Iexport=0
      DO ifield=1,Nexport(Iocean)
        id=ExportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)

        SELECT CASE (TRIM(code))

          CASE ('Pair')                   ! surface air pressure

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%PSFC(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Hair')                   ! surface air humidity

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
!
!  Calculate pressure at 2 m using hypsometric equation.
!  Assume temperature at 2 m as surface (0 m) temperature.
!
                cff1 = grid%PSFC(i,j) /                                 &
     &                 (exp((9.81*2.0)/(287.0*grid%T2(i,j))))
!
!  Compute specific humidity using the 2-m mixing ratio and 2-m pressure.
!
                rnum = grid%Q2(i,j)*cff1
                rden = (grid%Q2(i,j)*(1.-0.622)+0.622)
                cff2 = rnum/rden 
!
!  Compute saturation specific humidity using Bolton equation 10.
!
                c04 = 17.67*(grid%T2(i,j)-273.15)
                c05 = (grid%T2(i,j)-273.15) + 243.5
                cff3 = 6.112*exp(c04/c05)
!
                AA(ij)=REAL(cff2/cff3,r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Tair')                   ! surface (2m) air temperature

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%T2(i,j)-273.15,r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Uwind')                  ! U-wind (10m) component

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%U10(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Vwind')                  ! V-wind (10m) component

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%V10(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('cloud')                  ! cloud fraction

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%CLDFRA(i,1,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('rain')                   ! precipitation

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(0.001*                                      &
     &                      (grid%RAINCV(i,j)+grid%RAINNCV(i,j)),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('LWrad')                  ! longwave radiation

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%GLW(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('SWrad')                  ! shortwave radiation

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%SWDOWN(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Lheat')                  ! latent heat flux

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%LH(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Sheat')                  ! sensible heat flux

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%HFX(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Ustr')                   ! surface U-wind stress (Pa)

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                cff1=1.0/(grid%em_alt(i,1,j)+eps)
                cff2=1.0/((grid%em_u_2(i,1,j)**2+                       &
     &                     grid%em_v_2(i,1,j)**2)**0.5+eps)
                ij=ij+1
                AA(ij)=REAL(cff1*cff2*                                  &
     &                      (grid%UST(i,j)**2)*(grid%em_u_2(i,1,j)),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Vstr')                   ! surface V-wind stress (Pa)

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                cff1=1.0/(grid%em_alt(i,1,j)+eps)
                cff2=1.0/((grid%em_u_2(i,1,j)**2+                       &
     &                     grid%em_v_2(i,1,j)**2)**0.5+eps)
                ij=ij+1
                AA(ij)=REAL(cff1*cff2*                                  &
     &                      (grid%UST(i,j)**2)*(grid%em_v_2(i,1,j))
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)

        END SELECT
      END DO
!
!  Send fields to ocean model.
!
      IF (Iexport.gt.0) THEN
        inpbuffer(2)=my_wtime(wtime)
        CALL MCT_SEND (atm2ocn_AV, WRFtoROMS, MyError)
        SendTime=SendTime+my_wtime(wtime)-inpbuffer(2)
        IF (MyError.ne.0) THEN
          IF (MyRank) THEN
            WRITE (6,10) 'ocean model, MyError = ', MyError
          END IF
          RETURN
        END IF
      END IF

# ifdef SWAN_COUPLING
!
!-----------------------------------------------------------------------
!  Export fields from atmosphere (WRF) to wave (SWAN) model.
!-----------------------------------------------------------------------
!
!  Schedule sending fields to the wave model.
!
      DO ifield=1,Nexport(Iwaves)
        id=ExportID(Iwaves)%val(ifield)
        code=ADJUSTL(Fields(id)%code)

        SELECT CASE (TRIM(code))

          CASE ('Uwind')                  ! U-wind (10m) component

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%U10(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

          CASE ('Vwind')                  ! V-wind (10m) component

            ij=0
            DO j=jds,jde
              DO i=ids,ide
                ij=ij+1
                AA(ij)=REAL(grid%V10(i,j),r8)
              END DO
            END DO
            CALL AttrVect_importRAttr (atm2ocn_AV, TRIM(code), AA, Asize)
            Iexport=Iexport+1

        END SELECT
      END DO
!
!  Send ocean fields to wave model.
!
      IF (Iexport.gt.0) THEN
        inpbuffer(2)=my_wtime(wtime)
        CALL MCT_SEND (atm2wav_AV, WRFtoSWAN, MyError)
        SendTime=SendTime+my_wtime(wtime)-inpbuffer(2)
        IF (MyError.ne.0) THEN
          IF (MyRank.eq.0) THEN
            WRITE (6,10) 'wave model, MyError = ', MyError
          END IF
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Import fields from ocean model (ROMS) to ocean atmosphere (WRF).
!  Currently, both atmosphere and ocean model grids are the same.
!  We need to revisit this logic to allow interpolation.
!-----------------------------------------------------------------------
!
!  Schedule receiving fields from ocean model.
!
      inpbuffer(1)=my_wtime(wtime)
      CALL MCT_Recv (ocn2atm_AV, WRFtoROMS, MyError)
      RecvTime=RecvTime+my_wtime(wtime)-inpbuffer(1)
      IF (MyError.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) 'ocean model, MyError = ', MyError
        END IF
        exit_flag=2
        RETURN
      END IF
!
!  Receive fields from atmosphere model.
!
      Iimport=0
      DO ifield=1,Nimport(Iocean)
        id=ImportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)

        SELECT CASE (TRIM(code))

          CASE ('SST')                    ! sea surface temperature

          CALL AttrVect_exportRAttr (ocn2atm_AV, TRIM(code), AA, Asize)
          ij=0
          DO j=jds,jde
            DO i=ids,ide
              ij=ij+1
              grid%sst(i,j)=(AA(ij)+273.15)*(1-grid%landmask(i,j))+     &
     &                       grid%sst(i,j)*grid%landmask(i,j)
            END DO
          END DO

        END SELECT
      END DO
!
!  Deallocate communication arrays.
!
      deallocate (AA)
!!    deallocate (points)
!
 10   FORMAT (' ATM2OCN_COUPLING - error while sending fields to ',     &
     &        a, i4)
 20   FORMAT (' ATM2OCN_COUPLING - error while receiving fields from ', &
     &        a, i4)
 30   FORMAT (6x,'ATM2OCN   - (', i2.2, ') imported and (', i2.2,       &
     &        ') exported fields,', t62, 't = ', a,/, 16x,              &
     &        '- SWAN coupling exchanges wait clock (s):',/, 19x,       &
     &        '(Recv= ', 1p,e14.8,0p, ' Send= ', 1p,e14.8,0p,')')

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
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (WRFtoROMS, MyError)
      CALL AttrVect_clean (atm2ocn_AV, MyError)
      CALL GlobalSegMap_clean (GSMapWRF, MyError)

      END SUBROUTINE finalize_atm2ocn_coupling
#endif
      END MODULE atm_coupler_mod
