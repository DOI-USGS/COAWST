#include "cppdefs.h"
      SUBROUTINE memory
!
!svn $Id: memory.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes and report estimates of dynamic memory and    !
!  automatic memory requirements for current application.              !
!                                                                      !
!  The dynamical memory is that associated with the ocean state arrays,!
!  and it is allocated at runtime, and it is persistent until the ROMS !
!  termination of the execution.                                       !
!                                                                      !
!  The automatic arrays appear in subroutines and functions for        !
!  temporary local computations. They are created on entry to the      !
!  subroutine for intermediate computations and disappear on exit.     !
!  The automatic arrays (meaning non-static) are either allocated on   !
!  heap or stack memory.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
#endif
      USE mod_netcdf, ONLY : Matts, Mdims, Mvars, NvarD, NvarA
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ng, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: Nlevels, Ntiles
!
      real(r8) :: Avalue, bytefac, megabytefac, size1d, size2d
      real(r8) :: sumAsize, sumBsize, sumDsize
      real(r8) :: totalAsize, totalBsize, totalDsize
!
      real(r8), parameter :: spv = 0.0_r8
!
      real(r8), allocatable ::  Asize(:,:)      ! automatic arrays
      real(r8), allocatable ::  Bsize(:,:)      ! automatic mpi-buffers
      real(r8), allocatable ::  Dsize(:,:)      ! dynamic array
      real(r8), allocatable ::  IOsize(:,:)     ! NetCDF I/O
!
!-----------------------------------------------------------------------
!  Report estimate of dynamic memory and automatic memory requirements.
!-----------------------------------------------------------------------
!
!  Allocate and initialize.
!
!$OMP MASTER
      Ntiles=MAXVAL(NtileI)*MAXVAL(NtileJ)-1
      IF (.not.allocated(Asize)) THEN
        allocate ( Asize(0:Ntiles,Ngrids) )
        Asize=spv
      END IF
      IF (.not.allocated(Bsize)) THEN
        allocate ( Bsize(0:Ntiles,Ngrids) )
        Bsize=spv
      END IF
      IF (.not.allocated(Dsize)) THEN
        allocate ( Dsize(0:Ntiles,Ngrids) )
        Dsize=spv
      END IF
      IF (.not.allocated(IOsize)) THEN
        allocate ( IOsize(0:Ntiles,Ngrids) )
        IOsize=spv
      END IF
!
!  Determine size floating-point arrays in bytes.  We could use the
!  Fortran 2008 standard function STORAGE_SIZE.  However since ROMS
!  is double-precision by default, we just set its value to 64 bits
!  or 8 bytes (1 byte = 8 bits).  The number of array elements is
!  multiplied by the megabytes factor.
!
      bytefac=REAL(KIND(bytefac),r8)    ! r8 kind in bytes
      megabytefac=bytefac*1.0E-6_r8     ! 1 Mb = 1.0E+6 bytes (SI units)
!
!  Add static memory requirements for processing NetCDF data.  The
!  variables are declared in "mod_netcdf".  Notice that a single
!  character has a size of eight bits (1 byte).
!
      Dmem(1)=Dmem(1)+REAL(Matts,r8)               ! att_kind
      Dmem(1)=Dmem(1)+2.0_r8*REAL(Mdims,r8)        ! dim_id,dim_size
      Dmem(1)=Dmem(1)+5.0_r8*REAL(Mvars,r8)        ! var_*
      Dmem(1)=Dmem(1)+REAL(NvarD*Mvars,r8)         ! var_dim
      Dmem(1)=Dmem(1)+2.0_r8*REAL(NvarD,r8)        ! var_Dids,var_Dsize
      Dmem(1)=Dmem(1)+2.0_r8*REAL(NvarA,r8)        ! var_Aint,var_Afloat
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Matts,r8)   ! att_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Mdims,r8)   ! dim_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Mvars,r8)   ! dim_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*NvarA,r8)   ! var_Aname
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*NvarD,r8)   ! var_Dname
      Dmem(1)=Dmem(1)+0.125_r8*REAL(1024*NvarA,r8) ! var_Achar
!
!  Estimate automatic memory requirements (megabytes) by looking at the
!  routines that use it most, like step2d, step3d_t, or NetCDF I/O.
!  (see memory.txt for more information).
!
      DO ng=1,Ngrids
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
#ifdef NESTING
          IminS=BOUNDS(ng)%Istr(tile)-4
          ImaxS=BOUNDS(ng)%Iend(tile)+3
          JminS=BOUNDS(ng)%Jstr(tile)-4
          JmaxS=BOUNDS(ng)%Jend(tile)+3
#else
          IminS=BOUNDS(ng)%Istr(tile)-3
          ImaxS=BOUNDS(ng)%Iend(tile)+3
          JminS=BOUNDS(ng)%Jstr(tile)-3
          JmaxS=BOUNDS(ng)%Jend(tile)+3
#endif
          size1d=REAL((ImaxS-IminS+1),r8)
          size2d=REAL((ImaxS-IminS+1)*(JmaxS-JminS+1),r8)
#ifdef SOLVE3D
          Asize(tile,ng)=megabytefac*                                   &
     &                   (4.0_r8*size1d*REAL(N(ng)+1,r8)+               &
     &                    7.0_r8*size2d+                                &
     &                    5.0_r8*size2d*REAL(N(ng),r8)+                 &
     &                    1.0_r8*size2d*REAL(N(ng)*NT(ng),r8))
# if !(defined PARALLEL_IO && defined DISTRIBUTE)
#  ifdef INLINE_2DIO
          Nlevels=1
#  else
          Nlevels=N(ng)+1
#  endif
          IOsize(tile,ng)=megabytefac*2.0_r8*                           &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2)*(Nlevels),r8)
# else
          IOsize(tile,ng)=megabytefac*                                  &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2)*(N(ng)+1),r8)
# endif
#else
          Asize(tile,ng)=megabytefac*                                   &
     &                   (38.0_r8*size2d)
          IOsize(tile,ng)=megabytefac*                                  &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2),r8)
#endif
        END DO
      END DO
!
!  Determine total maximum value of dynamic-memory and automatic-memory
!  requirements, and convert number of array elements to megabytes.
!
#ifdef DISTRIBUTE
      Bsize(0:numthreads-1,1:Ngrids)=spv
      Dsize(0:numthreads-1,1:Ngrids)=spv
      DO ng=1,Ngrids
        Bsize(MyRank,ng)=BmemMax(ng)*1.0E-6_r8        ! already in bytes
        Dsize(MyRank,ng)=megabytefac*Dmem(ng)
        CALL mp_collect (ng, iNLM, numthreads, spv, Bsize(0:,ng))
        CALL mp_collect (ng, iNLM, numthreads, spv, Dsize(0:,ng))
      END DO
#else
      Dsize(0:numthreads-1,1:Ngrids)=spv
      DO ng=1,Ngrids
        Dsize(0,ng)=megabytefac*Dmem(ng)
      END DO
#endif
!
!  Report dynamic and automatic memory requirements.
!
      IF (Master) THEN
        WRITE (stdout,"(/,80('>'))")
        totalAsize=0.0_r8
        totalBsize=0.0_r8
        totalDsize=0.0_r8
        DO ng=1,Ngrids
          sumAsize=0.0_r8
          sumBsize=0.0_r8
          sumDsize=0.0_r8
#ifdef SOLVE3D
          WRITE (stdout,10) ng, Lm(ng), Mm(ng), N(ng),                  &
     &                      NtileI(ng), NtileJ(ng)
#else
          WRITE (stdout,10) ng, Lm(ng), Mm(ng),                         &
     &                      NtileI(ng), NtileJ(ng)
#endif
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            Avalue=MAX(Asize(tile,ng), Bsize(tile,ng), IOsize(tile,ng))
            sumAsize=sumAsize+Avalue
            sumBsize=sumBsize+Bsize(tile,ng)
            sumDsize=sumDsize+Dsize(tile,ng)
            WRITE (stdout,20) tile, Dsize(tile,ng), Avalue,             &
#ifdef DISTRIBUTE
     &                        Dsize(tile,ng)+Avalue, Bsize(tile,ng)
#else
     &                        Dsize(tile,ng)+Avalue
#endif
          END DO
          totalAsize=totalAsize+sumAsize
          totalBsize=totalBsize+sumBsize
          totalDsize=totalDsize+sumDsize
          IF (Ngrids.gt.1) THEN
            WRITE (stdout,30) '  SUM', sumDsize, sumAsize,              &
#ifdef DISTRIBUTE
     &                                 sumAsize+sumDsize, sumBsize
#else
     &                                 sumAsize+sumDsize
#endif
          ELSE
            WRITE (stdout,30) 'TOTAL', sumDsize, sumAsize,              &
#ifdef DISTRIBUTE
     &                                 sumAsize+sumDsize, sumBsize
#else
     &                                 sumAsize+sumDsize
#endif
          END IF
        END DO
        IF (Ngrids.gt.1) THEN
            WRITE (stdout,30) 'TOTAL', totalDsize, totalAsize,          &
#ifdef DISTRIBUTE
     &                                 totalAsize+totalDsize, totalBsize
#else
     &                                 totalAsize+totalDsize
#endif
        END IF
        WRITE (stdout,"(/,80('<'))")
      END IF
!
!  Deallocate dynamic and automatic memory local arrays.
!
      IF (allocated(Asize))  deallocate ( Asize )
      IF (allocated(Bsize))  deallocate ( Bsize )
      IF (allocated(Dsize))  deallocate ( Dsize )
      IF (allocated(IOsize)) deallocate ( IOsize )
!$OMP END MASTER
!$OMP BARRIER
!
#ifdef SOLVE3D
 10   FORMAT (/,' Dynamic and Automatic memory (MB) usage for Grid ',   &
     &          i2.2,':',2x,i0,'x',i0,'x',i0,2x,'tiling: ',i0,'x',i0,   &
     &          /,/,5x,'tile',10x,'Dynamic',8x,'Automatic',             &
     &          12x,'USAGE',                                            &
# ifdef DISTRIBUTE
     &          6x,'MPI-Buffers',                                       &
# endif
     &          /)
#else
 10   FORMAT (/,' Dynamic and Automatic memory (MB) usage for Grid ',   &
     &          i2.2,':',2x,i0,'x',i0,2x,'tiling: ',i0,'x',i0,          &
     &          /,/,5x,'tile',10x,'Dynamic',8x,'Automatic',             &
     &          12x,'USAGE',                                            &
# ifdef DISTRIBUTE
     &          6x,'MPI-Buffers',                                       &
# endif
     &          /)
#endif
 20   FORMAT (4x,i5,4(4x,f13.2))
 30   FORMAT (/,4x,a,4(4x,f13.2))
!
      RETURN
      END SUBROUTINE memory
