#include "cppdefs.h"
#ifdef FULL_GRID
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#else
# define IR_RANGE Istr,Iend
# define IU_RANGE IstrU,Iend
# define JR_RANGE Jstr,Jend
# define JV_RANGE JstrV,Jend
#endif

      MODULE wpoints_mod

#if defined PROPAGATOR || \
   (defined MASKING && (defined READ_WATER || defined WRITE_WATER))
!
!svn $Id: wpoints.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets variables associate with reading and writing of   !
!  water points data and/or determining size of packed state vector.   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: wpoints

      CONTAINS
!
!***********************************************************************
      SUBROUTINE wpoints (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL wpoints_tile (ng, tile, model,                               &
     &                   LBi, UBi, LBj, UBj,                            &
# ifdef MASKING
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  ifdef PROPAGATOR
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
#  endif
# endif
     &                   GRID(ng) % h)

      RETURN
      END SUBROUTINE wpoints
!
!***********************************************************************
      SUBROUTINE wpoints_tile (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj,                      &
# ifdef MASKING
     &                         pmask, rmask, umask, vmask,              &
#  ifdef PROPAGATOR
     &                         IJwaterR, IJwaterU, IJwaterV,            &
#  endif
# endif
     &                         h)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
# ifdef PROPAGATOR
      USE mod_storage
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_gather, mp_reduce
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
#   ifdef PROPAGATOR
      integer, intent(out) :: IJwaterR(LBi:,LBj:)
      integer, intent(out) :: IJwaterU(LBi:,LBj:)
      integer, intent(out) :: IJwaterV(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: h(LBi:,LBj:)
# else
#  ifdef MASKING
#   ifdef PROPAGATOR
      integer, intent(out) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(out) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(out) :: IJwaterV(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: my_Nxyp = 0
      integer :: my_Nxyr = 0
      integer :: my_Nxyu = 0
      integer :: my_Nxyv = 0
# ifdef PROPAGATOR
      integer :: my_NwaterR = 0
      integer :: my_NwaterU = 0
      integer :: my_NwaterV = 0
# endif

# ifdef PROPAGATOR
      integer :: Imin, Imax, Jmin, Jmax
# endif
      integer :: NSUB, Npts, i, ic, ij, j
# ifdef DISTRIBUTE
#  ifdef PROPAGATOR
      integer :: block_size

      real(r8), dimension(3) :: wp_buffer
      character (len=3), dimension(3) :: wp_handle
#  endif
#  if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
      real(r8), dimension(4) :: io_buffer
      character (len=3), dimension(4) :: io_handle
#  endif
# endif
# if defined READ_WATER || defined PROPAGATOR
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: mask
# endif

# include "set_bounds.h"

# if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
!
!-----------------------------------------------------------------------
!  Determine number of water points. 
!-----------------------------------------------------------------------
!
!  Determine interior and boundary points for IO purposes.
!
      my_Nxyr=0
      my_Nxyp=0
      my_Nxyu=0
      my_Nxyv=0
      DO j=MAX(0,JstrR),MIN(JendR,Mm(ng)+1)
        DO i=MAX(0,IstrR),MIN(IendR,Lm(ng)+1)
          IF (rmask(i,j).gt.0.0_r8) THEN
            my_Nxyr=my_Nxyr+1
          END IF
        END DO
        DO i=Istr,MIN(IendR,Lm(ng)+1)
          IF (umask(i,j).gt.0.0_r8) THEN
            my_Nxyu=my_Nxyu+1
          END IF
        END DO
      END DO
      DO j=Jstr,MIN(JendR,Mm(ng)+1) 
        DO i=Istr,MIN(IendR,Lm(ng)+1)
          IF (pmask(i,j).eq.1.0_r8) THEN
            my_Nxyp=my_Nxyp+1
          END IF
        END DO
        DO i=MAX(0,IstrR),MIN(IendR,Lm(ng)+1)
          IF (vmask(i,j).gt.0.0_r8) THEN
            my_Nxyv=my_Nxyv+1
          END IF
        END DO
      END DO
# endif
# ifdef PROPAGATOR
!
!-----------------------------------------------------------------------
!  Determine number of water points for propagator state vector.
!-----------------------------------------------------------------------
!
      my_NwaterR=0
      my_NwaterU=0
      my_NwaterV=0
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            my_NwaterR=my_NwaterR+1
          END IF
#  else
          my_NwaterR=my_NwaterR+1
#  endif
        END DO
        DO i=IU_RANGE
#  ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            my_NwaterU=my_NwaterU+1
          END IF
#  else
          my_NwaterU=my_NwaterU+1
#  endif
        END DO
      END DO
      DO j=JV_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            my_NwaterV=my_NwaterV+1
          END IF
#  else
          my_NwaterV=my_NwaterV+1
#  endif
        END DO
      END DO
# endif
!
!  Determine global number of points (parallel reduction).
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (MAX_WATER)
      IF (block_count.eq.0) THEN
# if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
        Nxyp(ng)=0
        Nxyr(ng)=0
        Nxyu(ng)=0
        Nxyv(ng)=0
# endif
# ifdef PROPAGATOR
        NwaterR(ng)=0
        NwaterU(ng)=0
        NwaterV(ng)=0
# endif
      END IF
# if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
      Nxyp(ng)=Nxyp(ng)+my_Nxyp
      Nxyr(ng)=Nxyr(ng)+my_Nxyr
      Nxyu(ng)=Nxyu(ng)+my_Nxyu
      Nxyv(ng)=Nxyv(ng)+my_Nxyv
# endif
# ifdef PROPAGATOR
      NwaterR(ng)=NwaterR(ng)+my_NwaterR
      NwaterU(ng)=NwaterU(ng)+my_NwaterU
      NwaterV(ng)=NwaterV(ng)+my_NwaterV
# endif
      block_count=block_count+1
      IF (block_count.eq.NSUB) THEN
        block_count=0
# ifdef DISTRIBUTE
#  if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
        io_buffer(1)=REAL(Nxyp(ng),r8)
        io_buffer(2)=REAL(Nxyr(ng),r8)
        io_buffer(3)=REAL(Nxyu(ng),r8)
        io_buffer(4)=REAL(Nxyv(ng),r8)
        io_handle(1)='SUM'
        io_handle(2)='SUM'
        io_handle(3)='SUM'
        io_handle(4)='SUM'
        CALL mp_reduce (ng, model, 4, io_buffer, io_handle)
        Nxyp(ng)=INT(io_buffer(1))
        Nxyr(ng)=INT(io_buffer(2))
        Nxyu(ng)=INT(io_buffer(3))
        Nxyv(ng)=INT(io_buffer(4))
#  endif
#  ifdef PROPAGATOR
        wp_buffer(1)=REAL(NwaterR(ng),r8)
        wp_buffer(2)=REAL(NwaterU(ng),r8)
        wp_buffer(3)=REAL(NwaterV(ng),r8)
        wp_handle(1)='SUM'
        wp_handle(2)='SUM'
        wp_handle(3)='SUM'
        CALL mp_reduce (ng, model, 3, wp_buffer, wp_handle)
        NwaterR(ng)=INT(wp_buffer(1))
        NwaterU(ng)=INT(wp_buffer(2))
        NwaterV(ng)=INT(wp_buffer(3))
#  endif
# endif
      END IF
!$OMP END CRITICAL (MAX_WATER)

# ifdef PROPAGATOR
!
!-----------------------------------------------------------------------
!  Determine size of state vector.
!-----------------------------------------------------------------------
#  if defined FORCING_SV
!
!  The state vector contains surface momentum stresses.
!
      Mstate(ng)=NwaterU(ng)+                                           &
     &           NwaterV(ng)
#  elif defined SO_SEMI
!
!  The state vector contains surface momentum stresses and surface
!  tracer fluxes.
!
      Mstate(ng)=0
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        Mstate(ng)=Mstate(ng)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        Mstate(ng)=Mstate(ng)+NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO i=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(i))) THEN
          Mstate(ng)=Mstate(ng)+NwaterR(ng)
        END IF
      END DO
#   endif
#  else
#   ifdef SOLVE3D
!
!  The state vector contains free-surface, 3D momentum components, and
!  tracers.
!
      Mstate(ng)=NwaterR(ng)+                                           &
     &           NwaterU(ng)*N(ng)+                                     &
     &           NwaterV(ng)*N(ng)+                                     &
     &           NwaterR(ng)*N(ng)*NT(ng)
#   else
!
!  The state vector contains free-surface and 2D momentum components.
!
      Mstate(ng)=NwaterR(ng)+                                           &
     &           NwaterU(ng)+                                           &
     &           NwaterV(ng)
#   endif
#  endif
#  ifdef DISTRIBUTE
!
!  Deternine state vector starting and ending blocking indices.
!
      Nstate(ng)=(Mstate(ng)+numnodes-1)/numnodes
      block_size=Nstate(ng)
      IF (Master) THEN
        WRITE (stdout,10) ng, Mstate(ng), Nstate(ng)
 10     FORMAT (' Size of FULL state arrays for propagator of grid ',   &
     &          i2.2,', Mstate = ', i10,/,38x,                          &
     &          'NODE partition, Nstate = ', i10,/,/,                   &
     &          9x,'Node',7x,'Nstr',7x,'Nend',7x,'Size',/)
        DO ic=0,numnodes-1
           i=1+ic*block_size
           j=MIN(Mstate(ng),i+block_size-1)
           WRITE (stdout,20) ic, i, j, j-i+1
 20        FORMAT (11x, i2, 3(1x,i10))
        END DO          
        WRITE (stdout,'(/)')
      END IF
      Nstr(ng)=1+MyRank*block_size
!!    Nend(ng)=MIN(Mstate(ng),Nstr(ng)+block_size-1)
      Nend(ng)=Nstr(ng)+block_size-1
      Mstate(ng)=Nstate(ng)*numnodes
#  else
      Nstr(ng)=1
      Nend(ng)=Mstate(ng)
      Nstate(ng)=Mstate(ng)
      IF (Master) THEN
        WRITE (stdout,10) ng, Mstate(ng)
 10     FORMAT (' Size of FULL state arrays for propagator of grid ',   &
     &          i2.2, ', Mstate = ', i10,/)
      END IF
#  endif
# endif

# if defined MASKING
!
!-----------------------------------------------------------------------
!  If distributed-memory, set process water indices arrays.
!-----------------------------------------------------------------------
#  ifdef DISTRIBUTE
!
!  Gather masking at RHO-points from all nodes.
!
      CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, 1, 1,              &
     &                r2dvar, 1.0_r8,                                   &
     &                rmask, rmask, Npts, mask)
      CALL mp_bcastf (ng, model, mask, Npts)
#  else
!
!  Load masking at RHO-points in mask vector.
!
      ic=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ic=ic+1
          mask(ic)=rmask(i,j)
        END DO
      END DO
      Npts=ic
#  endif
#  ifdef READ_WATER
!
!  Determine water RHO-points for IO purposes.
!
      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,r2dvar)=ij
        END IF
      END DO
      my_Nxyr=ic
#  endif
#  ifdef PROPAGATOR
!
!  Determine water points for the propagator state vector at RHO-points.
!
      ic=0      
      ij=0
#   ifdef FULL_GRID
      Imin=0
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
#   else
      Imin=1
      Imax=Lm(ng)
      Jmin=1
      Jmax=Mm(ng)
#   endif
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ij=ij+1
          IF ((mask(ij).gt.0.0_r8).and.                                 &
     &        (Imin.le.i).and.(i.le.Imax).and.                          &
     &        (Jmin.le.j).and.(j.le.Jmax)) THEN
            ic=ic+1
            mask(ij)=REAL(ic,r8)
          ELSE
            mask(ij)=0.0_r8
          END IF
        END DO
      END DO
!
      ij=0
      DO j=0,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ij=ij+1
          IF ((rILB(ng).le.i).and.(i.le.rIUB(ng)).and.                  &
     &        (rJLB(ng).le.j).and.(j.le.rJUB(ng))) THEN
            IF (mask(ij).gt.0.0_r8) THEN
              IJwaterR(i,j)=INT(mask(ij))
            ELSE
              IJwaterR(i,j)=0
            END IF
          END IF
        END DO
      END DO
#  endif

#  ifdef READ_WATER
#   ifdef DISTRIBUTE
!
!  Gather masking at PSI-points from all nodes.
!
      CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, 1, 1,              &
     &                p2dvar, 1.0_r8,                                   &
     &                pmask, pmask, Npts, mask)
      CALL mp_bcastf (ng, model, mask, Npts)
#   else
!
!  Load masking at PSI-points in mask vector.
!
      ic=0
      DO j=1,Mm(ng)
        DO i=1,Lm(ng)
          ic=ic+1
          mask(ic)=pmask(i,j)
        END DO
      END DO
      Npts=ic
#   endif
!
!  Determine water PSI-points for IO purposes.
!
      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,p2dvar)=ij
        END IF
      END DO
      my_Nxyp=ic
#  endif
#  ifdef DISTRIBUTE
!
!  Gather masking at U-points from all nodes.
!
      CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, 1, 1,              &
     &                u2dvar, 1.0_r8,                                   &
     &                umask, umask, Npts, mask)
      CALL mp_bcastf (ng, model, mask, Npts)
#  else
!
!  Load masking at U-points in mask vector.
!
      ic=0
      DO j=0,Mm(ng)+1
        DO i=1,Lm(ng)+1
          ic=ic+1
          mask(ic)=umask(i,j)
        END DO
      END DO
      Npts=ic
#  endif
#  ifdef READ_WATER
!
!  Determine water U-points for IO purposes.
!
      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,u2dvar)=ij
        END IF
      END DO
      my_Nxyu=ic
#  endif
#  ifdef PROPAGATOR
!
!  Determine water points for the propagator state vector at U-points.
!
      ic=0
      ij=0
#   ifdef FULL_GRID
      Imin=1
      Imax=Lm(ng)+1
      Jmin=0
      Jmax=Mm(ng)+1
#   else
      Imin=2
      Imax=Lm(ng)
      Jmin=1
      Jmax=Mm(ng)
#   endif
      DO j=0,Mm(ng)+1
        DO i=1,Lm(ng)+1
          ij=ij+1
          IF ((mask(ij).gt.0.0_r8).and.                                 &
     &        (Imin.le.i).and.(i.le.Imax).and.                          &
     &        (Jmin.le.j).and.(j.le.Jmax)) THEN
            ic=ic+1
            mask(ij)=REAL(ic,r8)
          ELSE
            mask(ij)=0.0_r8
          END IF
        END DO
      END DO
!
      ij=0
      DO j=0,Mm(ng)+1
        DO i=1,Lm(ng)+1
          ij=ij+1
          IF ((uILB(ng).le.i).and.(i.le.uIUB(ng)).and.                  &
     &        (uJLB(ng).le.j).and.(j.le.uJUB(ng))) THEN
            IF (mask(ij).gt.0.0_r8) THEN
              IJwaterU(i,j)=INT(mask(ij))
            ELSE
              IJwaterU(i,j)=0
            END IF
          END IF
        END DO
      END DO
#  endif
#  ifdef DISTRIBUTE
!
!  Gather masking at V-points from all nodes.
!
      CALL mp_gather (ng, model, LBi, UBi, LBj, UBj, 1, 1,              &
     &                v2dvar, 1.0_r8,                                   &
     &                vmask, vmask, Npts, mask)
      CALL mp_bcastf (ng, model, mask, Npts)
#  else
!
!  Load masking at V-points in mask vector.
!
      ic=0
      DO j=1,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ic=ic+1
          mask(ic)=vmask(i,j)
        END DO
      END DO
      Npts=ic
#  endif
#  ifdef READ_WATER
!
!  Determine water V-points for IO purposes.
!
      ic=0
      DO ij=1,Npts
        IF (mask(ij).gt.0.0_r8) THEN
          ic=ic+1
          SCALARS(ng)%IJwater(ic,v2dvar)=ij
        END IF
      END DO
      my_Nxyv=ic
#  endif
#  ifdef PROPAGATOR
!
!  Determine water points for the propagator state vector at V-points.
!
      ic=0
      ij=0
#   ifdef FULL_GRID
      Imin=0
      Imax=Lm(ng)+1
      Jmin=1
      Jmax=Mm(ng)+1
#   else
      Imin=1
      Imax=Lm(ng)
      Jmin=2
      Jmax=Mm(ng)
#   endif
      DO j=1,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ij=ij+1
          IF ((mask(ij).gt.0.0_r8).and.                                 &
     &        (Imin.le.i).and.(i.le.Imax).and.                          &
     &        (Jmin.le.j).and.(j.le.Jmax)) THEN
            ic=ic+1
            mask(ij)=REAL(ic,r8)
          ELSE
            mask(ij)=0.0_r8
          END IF
        END DO
      END DO
!
      ij=0
      DO j=1,Mm(ng)+1
        DO i=0,Lm(ng)+1
          ij=ij+1
          IF ((vILB(ng).le.i).and.(i.le.vIUB(ng)).and.                  &
     &        (vJLB(ng).le.j).and.(j.le.vJUB(ng))) THEN
            IF (mask(ij).gt.0.0_r8) THEN
              IJwaterV(i,j)=INT(mask(ij))
            ELSE
              IJwaterV(i,j)=0
            END IF
          END IF
        END DO
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE wpoints_tile
#endif
      END MODULE wpoints_mod
