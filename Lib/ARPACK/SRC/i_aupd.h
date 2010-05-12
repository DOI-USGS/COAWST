c
c     %-----------------------------------------------------%
c     | Private variables used by  _AUPD routines are saved |
c     ! in a common block to allow checkpointing. All these |
c     ! variables need to be saved and recovered during the |
c     ! checkpointing restart.                              |
c     %-----------------------------------------------------%
c
      integer
     &           bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritz, ritzi, ritzr
      common /i_aupd/
     &           bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritz, ritzi, ritzr
