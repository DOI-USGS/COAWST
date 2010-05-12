#include "cppdefs.h"
#if defined WEAK_CONSTRAINT && defined OBS_IMPACT
      SUBROUTINE rep_matrix (ng, model, outLoop, NinnLoop)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine estimates the inverse of the stabilized representer    !
!  matrix using the Lanczos vectors from the inner-loops.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
# endif
!
      implicit none
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, model, outLoop, NinnLoop
!
!  Local variable declarations.
!
      integer :: i, j, iobs, ivec, innLoop

      real(r8) :: zbet

      real(r8), dimension(NinnLoop) :: zu, zgam, zt
!
!-----------------------------------------------------------------------
!  Compute observation impact on the weak constraint 4DVAR data
!  assimilation systems.
!-----------------------------------------------------------------------
!
      MASTER_THREAD : IF (Master) THEN

        DO iobs=1,Ndatum(ng)
          ad_ObsVal(iobs)=0.0_r8
        END DO
        DO innLoop=1,NinnLoop
          zt(innLoop)=0.0_r8
        END DO
!
!  Multiply TLmodVal by the tranpose matrix of Lanczos vectors.
!  Note that the factor of 1/SQRT(ObsErr) is added to convert to
!  x-space.
!
        DO innLoop=1,NinnLoop
          DO iobs=1,Ndatum(ng)
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              zt(innLoop)=zt(innLoop)+ObsScale(iobs)*                   &
     &                    zcglwk(iobs,innLoop,outLoop)*TLmodVal(iobs)/  &
     &                    SQRT(ObsErr(iobs))
            END IF
          END DO
        END DO
!
!  Now multiply the result by the inverse tridiagonal matrix.
!
        zbet=cg_delta(1,outLoop)
        zu(1)=zt(1)/zbet
        DO ivec=2,NinnLoop
          zgam(ivec)=cg_beta(ivec,outLoop)/zbet
          zbet=cg_delta(ivec,outLoop)-cg_beta(ivec,outLoop)*zgam(ivec)
          zu(ivec)=(zt(ivec)-cg_beta(ivec,outLoop)*                     &
     &              zu(ivec-1))/zbet
        END DO

        DO ivec=NinnLoop-1,1,-1
          zu(ivec)=zu(ivec)-zgam(ivec+1)*zu(ivec+1)
        END DO
!
!  Finally multiply by the matrix of Lanczos vactors.
!  Note that the factor of 1/SQRT(ObsErr) is added to covert to
!  x-space.
!
        DO iobs=1,Ndatum(ng)
          DO innLoop=1,NinnLoop
            IF (ObsErr(iobs).NE.0.0_r8) THEN
              ad_ObsVal(iobs)=ad_ObsVal(iobs)+                          &
     &                        ObsScale(iobs)*                           &
     &                        zcglwk(iobs,innLoop,outLoop)*zu(innLoop)/ &
     &                        SQRT(ObsErr(iobs))
            END IF
          END DO
        END DO

      END IF MASTER_THREAD

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast new solution to other nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastf (ng, model, ad_ObsVal)
# endif
!
      RETURN
      END SUBROUTINE rep_matrix
#else
      SUBROUTINE rep_matrix
      RETURN
      END SUBROUTINE rep_matrix
#endif
