c
c     %---------------------------------------------------%
c     | Private variables used by  _AITR single precision |
c     | routines are saved in common blocks to facilitate |
c     | checkpointing.  All  these  variables  need to be |
c     | saved and recovered during checkpointing restart. |
c     %---------------------------------------------------%
c
      logical
     &           orth1, orth2, rstart, step3, step4
      common /lsaitr/
     &           orth1, orth2, rstart, step3, step4

      integer
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl
      common /isaitr/
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl

      real
     &           betaj, ovfl, rnorm1, safmin, smlnum, ulp, unfl, wnorm
      common /rsaitr/
     &           betaj, ovfl, rnorm1, safmin, smlnum, ulp, unfl, wnorm
