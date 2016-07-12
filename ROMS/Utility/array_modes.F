#include "cppdefs.h"
      MODULE array_modes_mod
#if defined WEAK_CONSTRAINT && (defined ARRAY_MODES || defined CLIPPING)
!
!svn $Id: array_modes.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used to process the requested eigenvector(s)     !
!  of the  stabilized representer matrix  when computing the array     !
!  modes or clipped spectrum analysis.                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Bennett, A.F, 2002: Inverse Modeling of the Ocean and Atmosphere,   !
!     Cambridge University Press, p 49-51.                             !
!                                                                      !
!======================================================================!
!
      implicit none

      PRIVATE
# ifdef ARRAY_MODES
      PUBLIC :: rep_check
      PUBLIC :: rep_eigen
# endif
# ifdef CLIPPING
      PUBLIC :: rep_clip
# endif

      CONTAINS

# ifdef ARRAY_MODES
      SUBROUTINE rep_check (ng, model, outLoop, NinnLoop)
!
!=======================================================================
!                                                                      !
!  This routine checks the dot-product of the array mode with the      !
!  corresponding eigenvector of the stabilized representer matrix.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#  endif
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: iobs, innLoop

      real(r8) :: zsum
!
!-----------------------------------------------------------------------
!  Check the dot-product of the array mode with the corresponding
!  eigenvector of the stabilized representer matrix.
!-----------------------------------------------------------------------
!
      MASTER_THREAD : IF (Master) THEN

        DO iobs=1,Ndatum(ng)
          ADModVal(iobs)=0.0_r8
        END DO
!
!  Multiply desired eigenvector of Lanczos tridiagonal matrix
!  by the Lanczos vectors to obtain the corresponding eigenvector
!  of the preconditioned stabilized representer matrix.
!
        DO iobs=1,Ndatum(ng)
          DO innLoop=1,NinnLoop
            ADmodVal(iobs)=ADmodVal(iobs)+                              &
     &                     cg_zv(innLoop,Nvct,outLoop)*                 &
     &                     zcglwk(iobs,innLoop,outLoop)
          END DO
        END DO
!
!  Now convert ADmodVal back to physical units.
!
        DO iobs=1,Ndatum(ng)
          IF (ObsErr(iobs).ne.0.0_r8) THEN
            ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
          END IF
        END DO
!
!  Now compute the dot-product of the final solution with the initial
!  eigenvector.
!
        zsum=0.0_r8
        DO iobs=1,Ndatum(ng)
          zsum=zsum+TLmodVal(iobs)*ADmodVal(iobs)
        END DO
!
!  Compare the dot-product with (cg_Ritz-1).
!
        WRITE (stdout,10) zsum, cg_Ritz(Nvct,outLoop)-1
 10     FORMAT (/,' REP CHECK: zsum = ', 1p, e14.7,2x,                  &
     &          'cg_Ritz-1 = ', 1p, e14.7)

      END IF MASTER_THREAD

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
#  endif

      RETURN
      END SUBROUTINE rep_check

      SUBROUTINE rep_eigen (ng, model, outLoop, NinnLoop)
!
!=======================================================================
!                                                                      !
!  This routine computes the specified eigenvector, Nvct, of the       !
!  stabilized representer matrix.                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#  endif
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: iobs, innLoop
!
!-----------------------------------------------------------------------
!  Compute specified eignevector of the stabilized representer matrix.
!-----------------------------------------------------------------------
!
      MASTER_THREAD : IF (Master) THEN

        DO iobs=1,Ndatum(ng)
          ADModVal(iobs)=0.0_r8
        END DO
!
!  Multiply desired eigenvector of Lanczos tridiagonal matrix
!  by the Lanczos vectors to obtain the corresponding eigenvector
!  of the preconditioned stabilized representer matrix.
!
        DO iobs=1,Ndatum(ng)
          DO innLoop=1,NinnLoop
            ADmodVal(iobs)=ADmodVal(iobs)+                              &
     &                     cg_zv(innLoop,Nvct,outLoop)*                 &
     &                     zcglwk(iobs,innLoop,outLoop)
          END DO
        END DO
!
!  Now convert ADmodVal back to physical units.
!
        DO iobs=1,Ndatum(ng)
          IF (ObsErr(iobs).ne.0.0_r8) THEN
            ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
          END IF
        END DO

      END IF MASTER_THREAD

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
#  endif

      RETURN
      END SUBROUTINE rep_eigen
# endif

# ifdef CLIPPING
      SUBROUTINE rep_clip (ng, model, outLoop, NinnLoop)
!
!=======================================================================
!                                                                      !
!  This routine performs clipping of the analysis by disgarding        !
!  potentially unphysical array modes.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#  endif
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: iobs, ivec, innLoop

      real(r8), dimension(NinnLoop) :: zu

      real(r8), dimension(Ndatum(ng)) :: innov, rsvec
!
!-----------------------------------------------------------------------
!  Clipp the analysis by disgarding potentially unphysical array modes.
!-----------------------------------------------------------------------
!
      MASTER_THREAD : IF (Master) THEN
!
!  First compute the dot-product of innovation vector with each
!  selected eigenvector of the stabilized representer matrix.
!  All eigenvectors < Nvct are disgarded.
!
        DO iobs=1,Ndatum(ng)
          innov(iobs)=ObsVal(iobs)-NLmodVal(iobs)
        END DO
        DO ivec=Nvct,Ninner
          DO iobs=1,Ndatum(ng)
            rsvec(iobs)=0.0_r8
          END DO
          DO iobs=1,Ndatum(ng)
            DO innLoop=1,NinnLoop
              rsvec(iobs)=rsvec(iobs)+                                  &
     &                    cg_zv(innLoop,ivec,outLoop)*                  &
     &                    zcglwk(iobs,innLoop,outLoop)
            END DO
          END DO
!
!  Convert RSVEC back to physical units.
!
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).ne.0.0_r8) THEN
              rsvec(iobs)=rsvec(iobs)/SQRT(ObsErr(iobs))
            END IF
          END DO
          zu(ivec)=0.0_r8
          DO iobs=1,Ndatum(ng)
            zu(ivec)=zu(ivec)+innov(iobs)*rsvec(iobs)
          END DO
        END DO
!
!  Second divide by the eigenvalues of the stabilized representer
!  matrix.
!
        DO ivec=Nvct,Ninner
          zu(ivec)=zu(ivec)/cg_Ritz(ivec,outLoop)
        END DO
!
!  Finally form the weighted sum of the selected eigenvectors of the
!  stabilized representer matrix.
!
        DO iobs=1,Ndatum(ng)
          ADModVal(iobs)=0.0_r8
        END DO
!
!  Multiply desired eigenvector of Lanczos tridiagonal matrix
!  by the Lanczos vectors to obtain the corresponding eigenvector
!  of the preconditioned stabilized representer matrix.
!
        DO ivec=Nvct,Ninner
          DO iobs=1,Ndatum(ng)
            DO innLoop=1,NinnLoop
              ADmodVal(iobs)=ADmodVal(iobs)+                            &
     &                       cg_zv(innLoop,ivec,outLoop)*               &
     &                       zcglwk(iobs,innLoop,outLoop)*              &
     &                       zu(ivec)
            END DO
          END DO
        END DO
!
!  Now convert ADmodVal back to physical units.
!
        DO iobs=1,Ndatum(ng)
          IF (ObsErr(iobs).ne.0.0_r8) THEN
            ADmodVal(iobs)=ADmodVal(iobs)/SQRT(ObsErr(iobs))
          END IF
        END DO

      END IF MASTER_THREAD

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ADmodVal)
#  endif

      RETURN
      END SUBROUTINE rep_clip
# endif
#endif
      END MODULE array_modes_mod
