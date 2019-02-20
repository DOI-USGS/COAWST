#include "cppdefs.h"
      MODULE waves_coupler_mod
#ifdef WAVES_OCEAN
!svn $Id: waves_coupler.F 817 2007-06-07 20:01:05Z jcwarner $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group      Hernan G. Arango   !
!   Licensed under a MIT/X style license                               !
!   See License_ROMS.txt                                               !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between SWAN   !
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
!      implicit none
!
      PRIVATE

      PUBLIC :: initialize_ocean_coupling
      PUBLIC :: ocean_coupling
      PUBLIC :: finalize_ocean_coupling

      include 'mpif.h'
!
!  Declarations.
!
      TYPE(GlobalSegMap) :: GSMapSWAN         ! GloabalSegMap variables
 
      TYPE(AttrVect) :: ToOceanAV             ! AttrVect variables
      TYPE(AttrVect) :: FromOceanAV
 
      type(Router) :: RoutSWAN                ! Router variables

      CONTAINS

      SUBROUTINE initialize_ocean_coupling
!
!=======================================================================
!                                                                      !
!  Initialize waves and ocean models coupling stream.  This is the     !
!  training phase use to constuct  MCT  parallel interpolators and     !
!  stablish communication patterns.                                    !
!                                                                      !
!=======================================================================
!
#include "param.h"
#include "common.h"
#include "pass.h"
!
      include 'mpif.h'
!
!  Local variable declarations.  
!
      integer :: MyColor, MyCOMM, MyError, MyRank, MyKey, MyValue
      integer :: npoints, Asize, nprocs, localsize
      integer :: j, Isize, Jsize

      integer, pointer :: start(:), length(:)

!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
! Get local rank and size.
!
      CALL mpi_comm_rank (WAV_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (WAV_COMM_WORLD, nprocs, MyError)
!
! Initialize MCTworld.
!
      CALL MCTWorld_init (ncomps,MPI_COMM_WORLD,WAV_COMM_WORLD,WavId)
!
! Initialize a Global Segment Map for non-haloed transfer of data out of
! SWAN. Determine non-haloed start and length arrays for this processor.
!
      IF (nprocs.eq.1) THEN
        Isize=mr
        Jsize=nr
      ELSE
        write(*,*) 'Refdif can not be tiled !!'
      END IF
!
      allocate ( start(Jsize) )
      allocate ( length(Jsize) )
!
      DO j=1,Jsize
        length(j)=Isize
        start(j)=(j-1)*Isize+1
      END DO
!
      Asize=Isize*Jsize
!  
      CALL GlobalSegMap_init(GSMapSWAN,start,length,0,                  &
     &                       WAV_COMM_WORLD,WavId)
!
! Initialize Attribute Vector ToOceanAv to hold the data sent to ROMS.
!
      CALL AttrVect_init (ToOceanAV,                                    &
     &     rList="DWAVE:HWAVE:LWAVE:WAVE_BREAK:WAVE_DISSIP",            &
     &     lsize=Asize)
      CALL AttrVect_zero (ToOceanAV)
!
! Initialize Attribute Vector FromOceanAV that will have ROMS data in
! it.
!
      CALL AttrVect_init (FromOceanAV,                                  &
     &                    rList="XR:YR:DEPTH:UBAR:VBAR:ZETA",           &
     &                    lsize=Asize)
      CALL AttrVect_zero (FromOceanAV)
!
! Initialize a router to the Waves component.
!
      CALL Router_init (OcnId,GSMapSWAN,WAV_COMM_WORLD,RoutSWAN)
!
      deallocate (start)
      deallocate (length)

      RETURN
      END SUBROUTINE initialize_ocean_coupling

      SUBROUTINE ocean_coupling
!
!=======================================================================
!                                                                      !
!  This subroutine reads and writes the coupling data streams between  !
!  ocean and wave models. Currently, the following data streams are    !
!  processed:                                                          !
!  Fields sent to the OCEAN Model:                                     !
!                                                                      !
!     * Dwave      Wave direction.                                     !
!     * Hwave      Wave height.                                        !
!     * Lwave      Wave length.                                        !
!     * Wave_break Percent of breakig waves.                           !
!     * Wave_dissip Wave energy dissipation.                           !
!                                                                      !
!  Fields acquired from the OCEAN Model:                               !
!                                                                      !
!     * xr         x-rho coordinates                                   !
!     * yr         y-rho coordinates                                   !
!     * h          Bottom elevation.                                   !
!     * ubar       Depth integrated xi-direction velocity.             !
!     * vbar       Depth integrated eta-direction velocity.            !
!     * zeta       Water surface elevation.                            !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!      implicit none
!
#include "param.h"
#include "common.h"
#include "pass.h"
!
!
!  Imported variable declarations.
!
!      integer :: IX, IY
!
!  Local variable declarations.
!
      integer :: MyStatus, i, j, Asize, ierr, MyRank
      integer :: MyError, MySize, indx, Istr, Iend, Jstr, Jend
      real :: cff

      integer, pointer :: points(:)
      real(r8), pointer :: AA(:)
!
!-----------------------------------------------------------------------
!  Send wave fields to ROMS.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK (WAV_COMM_WORLD, MyRank, MyError)
      CALL MPI_COMM_SIZE (WAV_COMM_WORLD, NPROCS, MyError)
!
!  Size is the number of grid point on this processor.
!
      Asize=GlobalSegMap_lsize(GSMapSWAN,WAV_COMM_WORLD)
!
!  Load wave data into Attribute Vector Array ToOceanAV.
!
      allocate (AA(Asize))
      allocate (points(Asize))
      avdata=0.0
      points=0
!
!  Ask for points in this tile.
!
      CALL GlobalSegMap_Ordpnts (GSMapSWAN,MyRank,points)
!
!-----------------------------------------------------------------------
!  Schedule and send required fields to ocean  model.
!-----------------------------------------------------------------------
!
!  Wave direction.
!
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          AA(ij)=pass_theta(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToOceanAV, "DWAVE", AA, Asize)
!
!  Wave height.
!
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          AA(ij)=pass_height(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToOceanAV, "HWAVE", AA, Asize)
!
!  Wave number.
!
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          AA(ij)=pass_wavenum(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToOceanAV, "LWAVE", AA, Asize)
!
!  Wave break area.
!
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          AA(ij)=pass_area(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToOceanAV, "WAVE_BREAK", AA, Asize)
!
!  Wave dissipation.
!
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          AA(ij)=pass_diss(i,j)
        END DO
      END DO
      CALL AttrVect_importRAttr (ToOceanAV, "WAVE_DISSIP", AA, Asize)
!
!-----------------------------------------------------------------------
!  Send wave parameters to ROMS.
!-----------------------------------------------------------------------
!
        CALL MCT_SEND (ToOceanAV,RoutSWAN,MyError)
        IF (MYRANK.EQ.0) THEN
          WRITE (*,*) '== SWAN sent wave fields and Myerror= ',MyError
        ENDIF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling send fail swancplr, MyStatus= ', MyError
          CALL finalize_ocean_coupling ("Coupling failed swancplr")
        END IF
!
!-----------------------------------------------------------------------
!  Receive from ROMS: Depth, Water Level, VELX, and VELY.
!-----------------------------------------------------------------------
!
        CALL MCT_Recv (FromOceanAV, RoutSWAN, MyError)
        IF (MYRANK.EQ.0) THEN
          WRITE (*,*) '== SWAN recvd ocean fields and Myerror??= ',     &
     &                MyError
        END IF
        IF (MyError.ne.0) THEN
          WRITE (*,*) 'coupling fail swancplr, MyStatus= ', MyError
          CALL finalize_ocean_coupling ("Coupling failed swancplr")
        END IF
!
!  x-coordinate.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "XR", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          x_wave(i,j)=AA(ij)
        END DO
      END DO
!
!  y-coordinate.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "YR", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          y_wave(i,j)=AA(ij)
        END DO
      END DO
!
!  depth.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "DEPTH", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          depth_wave(i,j)=AA(ij)
        END DO
      END DO
!
!  ubar.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "UBAR", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          intp_U_wave(i,j)=AA(ij)
        END DO
      END DO
!
!  vbar.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "VBAR", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          intp_V_wave(i,j)=AA(ij)
        END DO
      END DO
!
!  zeta.
!
      CALL AttrVect_exportRAttr (FromOceanAV, "ZETA", AA, Asize)
      ij=0
      DO j=1,NR
        DO i=1,MR
          ij=ij+1
          intp_eta_wave(i,j)=AA(ij)
        END DO
      END DO
!
      deallocate (AA, points)

      RETURN
      END SUBROUTINE ocean_coupling

      SUBROUTINE finalize_ocean_coupling (string)
!
!=======================================================================
!                                                                    ===
!  This routines terminates execution during coupling error.         ===
!                                                                    ===
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: string
!
!  Local variable declarations.
!
      integer :: MyStatus
!
!-----------------------------------------------------------------------
!  Terminate MPI execution environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (RoutSWAN)
      CALL AttrVect_clean (ToOceanAV)
      CALL AttrVect_clean (FromOceanAV)
      CALL GlobalSegMap_clean (GSMapSWAN)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyStatus)

      STOP
      END SUBROUTINE finalize_ocean_coupling
#endif
      END MODULE waves_coupler_mod
