c
c     %---------------------------------------------------%
c     | Private variables used by  _AUP2 double precision |
c     | routines are saved in common blocks to facilitate |
c     | checkpointing.  All  these  variables  need to be |
c     | saved and recovered during checkpointing restart. |
c     %---------------------------------------------------%
c
      logical
     &           cnorm, getv0, initv, update, ushift
      common /ldaup2/
     &           cnorm, getv0, initv, update, ushift

      integer
     &           iter, kplusp, msglvl, nconv, nevbef, nev0, np0, numcnv
      common /idaup2/
     &           iter, kplusp, msglvl, nconv, nevbef, nev0, np0, numcnv

      double precision
     &           rnorm, eps23
      common /rdaup2/
     &           rnorm, eps23
