#include "cppdefs.h"
      MODULE math_tools
#ifdef INWAVE_SWAN_COUPLING
!
!svn $Id: math_tools.F 830 2017-01-24 21:21:11Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!    Submitted by John C. Warner and Maitane Olabarrieta               !
!=======================================================================
!                                                                      !
!  MODULE CONTAINS:                                                    !
!  This routine perfoms the Hilbert transform and some matrix          !
!  manipulatoins.                                                      !
!     1) Hilbert transform, calls fftpack                              !
!=======================================================================
!
      IMPLICIT NONE

      PRIVATE
      PUBLIC:: fftkind, hilbert

      INTEGER, PARAMETER:: fftkind = KIND(0.0d0)

      CONTAINS
!
!***********************************************************************
      subroutine hilbert(xi,m)
!***********************************************************************
!
!     Calls fftpack.
!     CFFT1I initializes the transform,
!     CFFT1F does a forward transform;
!     CFFT1B does a backward transform.
!
      IMPLICIT NONE

      integer                             :: m, n, ier, lensav, lenwrk
      complex(fftkind),dimension(m)                :: xi
      complex(fftkind),dimension(:),allocatable    :: x
      integer,dimension(:),allocatable    :: h
      real                                :: p
      real(kind=8),dimension(:),allocatable    :: work, wsave

      p=log(real(m))/log(2.)
      n=ceiling(p)
      n=2**n

      allocate(x(n))
      x=0
      x(1:m)=real(xi)

      lenwrk=2*n
      lensav=2*n+int(log(real(n,kind=8))/log(2.0D+00))+4
      allocate ( work(1:lenwrk) )
      allocate ( wsave(1:lensav) )
      work=0.
      wsave=0.

      call cfft1i ( n, wsave, lensav, ier )
      call cfft1f ( n, 1, x, n, wsave, lensav, work, lenwrk, ier )

      x=x*sqrt(real(n))           ! Scale factor

      allocate(h(n))

      h(1)=1
      h(2:(n/2))=2
      h((n/2)+1)=1
      h((n/2)+2:n)=0

      x=x*h

      deallocate(h)

      call cfft1b ( n, 1, x, n, wsave, lensav, work, lenwrk, ier )
      x=x/sqrt(real(n))           ! Scale factor

      xi=x(1:m)
      deallocate(x)
      deallocate(work, wsave)

      return

      end subroutine hilbert

#endif
      END MODULE math_tools
