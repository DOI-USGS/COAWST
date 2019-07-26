c
c     %---------------------------------------------------%
c     | Private variables used by  _AUP2 single precision |
c     | routines are saved in common blocks to facilitate |
c     | checkpointing.  All  these  variables  need to be |
c     | saved and recovered during checkpointing restart. |
c     %---------------------------------------------------%
c
      logical
     &           cnorm, getv0, initv, update, ushift
      common /lsaup2/
     &           cnorm, getv0, initv, update, ushift

      integer
     &           iter, kplusp, msglvl, nconv, nevbef, nev0, np0, numcnv
      common /isaup2/
     &           iter, kplusp, msglvl, nconv, nevbef, nev0, np0, numcnv

      real
     &           rnorm, eps23
      common /rsaup2/
     &           rnorm, eps23
