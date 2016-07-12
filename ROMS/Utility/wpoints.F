#include "cppdefs.h"
#ifdef FULL_GRID
# define IR_RANGE IstrT,IendT
# define IU_RANGE IstrP,IendT
# define JR_RANGE JstrT,JendT
# define JV_RANGE JstrP,JendT
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
!svn $Id: wpoints.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets variables associated with reading and writing of  !
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
     &                   IminS, ImaxS, JminS, JmaxS,                    &
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
     &                         IminS, ImaxS, JminS, JmaxS,              &
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
      USE distribute_mod, ONLY : mp_bcastf, mp_gather2d, mp_reduce
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
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
      integer :: my_Nxyp, my_Nxyr, my_Nxyu, my_Nxyv
# ifdef PROPAGATOR
      integer :: my_NwaterR, my_NwaterU, my_NwaterV
# endif

# ifdef PROPAGATOR
      integer :: Imin, Imax, Jmin, Jmax
# endif
      integer :: Uoff, Voff
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
!  Determine interior points for IO purposes.
!
      my_Nxyr=0
      my_Nxyp=0
      my_Nxyu=0
      my_Nxyv=0

      DO j=Jstr,Jend
        DO i=Istr,Iend
          IF (pmask(i,j).gt.0.0_r8) my_Nxyp=my_Nxyp+1
          IF (rmask(i,j).gt.0.0_r8) my_Nxyr=my_Nxyr+1
          IF (umask(i,j).gt.0.0_r8) my_Nxyu=my_Nxyu+1
          IF (vmask(i,j).gt.0.0_r8) my_Nxyv=my_Nxyv+1
        END DO
      END DO
!
!  Determine boundary points for IO purposes.  The assigment
!  below is done to account for periodicity.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        DO j=Jstr,Jend
          IF (rmask(Istr-1,j).gt.0.0_r8) my_Nxyr=my_Nxyr+1
          IF (vmask(Istr-1,j).gt.0.0_r8) my_Nxyv=my_Nxyv+1
        END DO
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        DO j=Jstr,Jend
          IF (pmask(Iend+1,j).gt.0.0_r8) my_Nxyp=my_Nxyp+1
          IF (rmask(Iend+1,j).gt.0.0_r8) my_Nxyr=my_Nxyr+1
          IF (umask(Iend+1,j).gt.0.0_r8) my_Nxyu=my_Nxyu+1
          IF (vmask(Iend+1,j).gt.0.0_r8) my_Nxyv=my_Nxyv+1
        END DO
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        DO i=Istr,Iend
          IF (rmask(i,Jstr-1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
          IF (umask(i,Jstr-1).gt.0.0_r8) my_Nxyu=my_Nxyu+1
        END DO
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        DO i=Istr,Iend
          IF (pmask(i,Jend+1).gt.0.0_r8) my_Nxyp=my_Nxyp+1
          IF (rmask(i,Jend+1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
          IF (umask(i,Jend+1).gt.0.0_r8) my_Nxyu=my_Nxyu+1
          IF (vmask(i,Jend+1).gt.0.0_r8) my_Nxyv=my_Nxyv+1
        END DO
      END IF
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        IF (rmask(Istr-1,Jstr-1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        IF (rmask(Iend+1,Jstr-1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
        IF (umask(Iend+1,Jstr-1).gt.0.0_r8) my_Nxyu=my_Nxyu+1
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        IF (rmask(Istr-1,Jend+1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
        IF (vmask(Istr-1,Jend+1).gt.0.0_r8) my_Nxyv=my_Nxyv+1
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        IF (pmask(Iend+1,Jend+1).gt.0.0_r8) my_Nxyp=my_Nxyp+1
        IF (rmask(Iend+1,Jend+1).gt.0.0_r8) my_Nxyr=my_Nxyr+1
        IF (umask(Iend+1,Jend+1).gt.0.0_r8) my_Nxyu=my_Nxyu+1
        IF (vmask(Iend+1,Jend+1).gt.0.0_r8) my_Nxyv=my_Nxyv+1
      END IF
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
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
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
# if defined MASKING && (defined READ_WATER || defined WRITE_WATER)
        IOBOUNDS(ng) % xy_psi = Nxyp(ng)
        IOBOUNDS(ng) % xy_rho = Nxyr(ng)
        IOBOUNDS(ng) % xy_u   = Nxyu(ng)
        IOBOUNDS(ng) % xy_v   = Nxyv(ng)
# endif
      END IF
!$OMP END CRITICAL (MAX_WATER)

# ifdef PROPAGATOR
!
!-----------------------------------------------------------------------
!  Determine size of state vector.
!-----------------------------------------------------------------------
#  if defined FORCING_SV || \
      defined SO_SEMI    || defined STOCHASTIC_OPT
!
!  The state vector may contain free-surface, momentum, tracers,
!  surface momentum stresses, amd surface tracer fluxes.
!
      Mstate(ng)=0
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        Mstate(ng)=Mstate(ng)+NwaterR(ng)
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        Mstate(ng)=Mstate(ng)+N(ng)*NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        Mstate(ng)=Mstate(ng)+N(ng)*NwaterV(ng)
      END IF
      DO i=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(i))) THEN
          Mstate(ng)=Mstate(ng)+N(ng)*NwaterR(ng)
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        Mstate(ng)=Mstate(ng)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        Mstate(ng)=Mstate(ng)+NwaterV(ng)
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        Mstate(ng)=Mstate(ng)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        Mstate(ng)=Mstate(ng)+NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO i=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(i))) THEN
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
      Nstate(ng)=(Mstate(ng)+numthreads-1)/numthreads
      block_size=Nstate(ng)
      IF (Master) THEN
        WRITE (stdout,10) ng, Mstate(ng), Nstate(ng)
 10     FORMAT (/,' Size of FULL state arrays for propagator of grid ', &
     &          i2.2,', Mstate = ', i10,/,38x,                          &
     &          'NODE partition, Nstate = ', i10,/,/,                   &
     &          9x,'Node',7x,'Nstr',7x,'Nend',7x,'Size',/)
        DO ic=0,numthreads-1
           i=1+ic*block_size
           j=MIN(Mstate(ng),i+block_size-1)
           WRITE (stdout,20) ic, i, j, j-i+1
 20        FORMAT (11x, i2, 3(1x,i10))
        END DO
        WRITE (stdout,'(/)')
      END IF
      Nstr(ng)=1+MyRank*block_size
      Nend(ng)=Nstr(ng)+block_size-1
      Mstate(ng)=Nstate(ng)*numthreads
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

# if ((defined READ_WATER && defined DISTRIBUTE) || \
      defined PROPAGATOR) && defined MASKING
!
!-----------------------------------------------------------------------
!  If distributed-memory, set process water indices arrays.
!-----------------------------------------------------------------------
!
!  Set offsets for U- and V-variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF

#  ifdef DISTRIBUTE
!
!  Gather masking at RHO-points from all nodes.
!
      CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                  &
     &                  0, r2dvar, 1.0_r8,                              &
     &                  rmask, rmask, Npts, mask, .FALSE.)
      CALL mp_bcastf (ng, model, mask)
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
#  if defined READ_WATER && defined DISTRIBUTE
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

#  if defined READ_WATER && defined DISTRIBUTE
!
!  Gather masking at PSI-points from all nodes.
!
      CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                  &
     &                  0, p2dvar, 1.0_r8,                              &
     &                  pmask, pmask, Npts, mask, .FALSE.)
      CALL mp_bcastf (ng, model, mask)
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
      CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                  &
     &                  0, u2dvar, 1.0_r8,                              &
     &                  umask, umask, Npts, mask, .FALSE.)
      CALL mp_bcastf (ng, model, mask)
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
#  if defined READ_WATER && defined DISTRIBUTE
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
      Imin=1+Uoff
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
      CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,                  &
     &                  0, v2dvar, 1.0_r8,                              &
     &                  vmask, vmask, Npts, mask, .FALSE.)
      CALL mp_bcastf (ng, model, mask)
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
#  if defined READ_WATER && defined DISTRIBUTE
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
      Jmin=1+Voff
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

# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE

      RETURN
      END SUBROUTINE wpoints_tile
#endif
      END MODULE wpoints_mod
