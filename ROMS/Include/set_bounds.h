/*
** svn $Id: set_bounds.h 814 2008-10-29 01:42:17Z jcwarner $
************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
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
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
#if defined COMPOSED_GRID || defined REFINED_GRID
      integer :: IstrC, IendC, JstrC, JendC
      integer :: IstrTU, JstrTV, IendTU, JendTV
#endif
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)

      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
#ifdef REFINED_GRID
      IstrTU=BOUNDS(ng)%IstrTU(tile)
      JstrTV=BOUNDS(ng)%JstrTV(tile)
      IendTU=BOUNDS(ng)%IendTU(tile)
      JendTV=BOUNDS(ng)%JendTV(tile)
#endif
