/*
** git $Id$
************************************************************************
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**
** The lower and upper bounds are allocated and assigned in "inp_par.F"
** by calling "get_tile" which is located in the "get_bounds.F" file.
**
*/

!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =xtr_BOUNDS(ng) % Istr   (tile)
      IstrB  =xtr_BOUNDS(ng) % IstrB  (tile)
      IstrM  =xtr_BOUNDS(ng) % IstrM  (tile)
      IstrP  =xtr_BOUNDS(ng) % IstrP  (tile)
      IstrR  =xtr_BOUNDS(ng) % IstrR  (tile)
      IstrT  =xtr_BOUNDS(ng) % IstrT  (tile)
      IstrU  =xtr_BOUNDS(ng) % IstrU  (tile)
      Iend   =xtr_BOUNDS(ng) % Iend   (tile)
      IendB  =xtr_BOUNDS(ng) % IendB  (tile)
      IendP  =xtr_BOUNDS(ng) % IendP  (tile)
      IendR  =xtr_BOUNDS(ng) % IendR  (tile)
      IendT  =xtr_BOUNDS(ng) % IendT  (tile)

      Jstr   =xtr_BOUNDS(ng) % Jstr   (tile)
      JstrB  =xtr_BOUNDS(ng) % JstrB  (tile)
      JstrM  =xtr_BOUNDS(ng) % JstrM  (tile)
      JstrP  =xtr_BOUNDS(ng) % JstrP  (tile)
      JstrR  =xtr_BOUNDS(ng) % JstrR  (tile)
      JstrT  =xtr_BOUNDS(ng) % JstrT  (tile)
      JstrV  =xtr_BOUNDS(ng) % JstrV  (tile)
      Jend   =xtr_BOUNDS(ng) % Jend   (tile)
      JendB  =xtr_BOUNDS(ng) % JendB  (tile)
      JendP  =xtr_BOUNDS(ng) % JendP  (tile)
      JendR  =xtr_BOUNDS(ng) % JendR  (tile)
      JendT  =xtr_BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =xtr_BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =xtr_BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =xtr_BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=xtr_BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=xtr_BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =xtr_BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =xtr_BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=xtr_BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =xtr_BOUNDS(ng) % Iendp3 (tile)            ! Iend+3

      Jstrm3 =xtr_BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =xtr_BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =xtr_BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=xtr_BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=xtr_BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =xtr_BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =xtr_BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=xtr_BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =xtr_BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
