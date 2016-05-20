subroutine SwanPunCollect ( blkndc )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2015  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!PUN!
!PUN!   Authors
!PUN!
!PUN!   41.36: Marcel Zijlema
!PUN!   41.40: Sander Hulst
!PUN!
!PUN!   Updates
!PUN!
!PUN!   41.36,     July 2012: New subroutine
!PUN!   41.40, November 2012: add global grid coordinates
!PUN!
!PUN!   Purpose
!PUN!
!PUN!   Determines for each grid vertex in global grid the node number
!PUN!   Also determines grid coordinates in global grid
!PUN!
!PUN!   Method
!PUN!
!PUN!   Collect specific data from all nodes
!PUN!
!PUN!   Modules used
!PUN!
!PUN    use mpi
!PUN    use ocpcomm4
!PUN    use swcomm2, only: XOFFS, YOFFS
!PUN    use m_parall, only: IAMMASTER
!PUN    use SwanGriddata, only: ivertg, nverts, nvertsg, xcugrd, ycugrd, xcugrdgl, ycugrdgl
!PUN    use SIZES, only: MYPROC, MNPROC
!PUN    use GLOBAL, only: COMM
!PUN!
!PUN    implicit none
!PUN!
!PUN!   Argument variables
!PUN!
!PUN    real,    dimension(nvertsg), intent(out) :: blkndc   ! array giving node number in each grid vertex in global grid
!PUN!
!PUN!   Local variables
!PUN!
!PUN    integer, save                            :: ient = 0 ! number of entries in this subroutine
!PUN    integer                                  :: ierr     ! error value of MPI call
!PUN    integer                                  :: j        ! loop counter
!PUN    integer                                  :: k        ! loop counter
!PUN    integer                                  :: nownv    ! number of vertices in own subdomain (without ghost vertices)
!PUN    !
!PUN    integer, dimension(:), allocatable       :: iarr     ! auxiliary integer array to gather data
!PUN    integer, dimension(:), allocatable       :: icount   ! array specifying array size of data received from each processor
!PUN    integer, dimension(:), allocatable       :: idsplc   ! array specifying the starting address of the incoming data from each processor, relative to the global array
!PUN    integer, dimension(:), allocatable       :: ivertp   ! vertex index of global grid in own subdomain (without ghost vertices)
!PUN    !
!PUN    real   , dimension(:), allocatable       :: arr      ! auxiliary real array to gather data
!PUN    real   , dimension(:), allocatable       :: blknd    ! node number per subdomain (without ghost vertices)
!PUN    real   , dimension(:), allocatable       :: xpl      ! user coordinates grid points (without ghost vertices)
!PUN    real   , dimension(:), allocatable       :: ypl      ! user coordinates grid points (without ghost vertices)
!PUN    !
!PUN    character(80)                            :: msgstr   ! string to pass message
!PUN!
!PUN!   Structure
!PUN!
!PUN!   Description of the pseudo code
!PUN!
!PUN!   Source text
!PUN!
!PUN    if (ltrace) call strace (ient,'SwanPunCollect')
!PUN    !
!PUN    nownv = count(ivertg>0)
!PUN    !
!PUN    allocate(ivertp(nownv))
!PUN    allocate( blknd(nownv))
!PUN    allocate(   xpl(nownv))
!PUN    allocate(   ypl(nownv))
!PUN    !
!PUN    ! determine node number per subdomain
!PUN    !
!PUN    k = 0
!PUN    do j = 1, nverts
!PUN       if ( ivertg(j) > 0 ) then
!PUN          k = k + 1
!PUN          ivertp(k) = ivertg(j)
!PUN          blknd (k) = real(MYPROC+1)
!PUN          xpl   (k) = xcugrd(j) + XOFFS
!PUN          ypl   (k) = ycugrd(j) + YOFFS
!PUN       endif
!PUN    enddo
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       allocate(icount(0:MNPROC-1))
!PUN       allocate(idsplc(0:MNPROC-1))
!PUN    endif
!PUN    !
!PUN    ! gather the array sizes to the master
!PUN    !
!PUN    call MPI_GATHER( nownv, 1, MPI_INTEGER, icount, 1, MPI_INTEGER, 0, COMM, ierr )
!PUN    if ( ierr /= MPI_SUCCESS ) then
!PUN       write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
!PUN       call msgerr( 4, trim(msgstr) )
!PUN       return
!PUN    endif
!PUN    !
!PUN    ! check consistency with respect to size of gathered data
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       if ( sum(icount) /= nvertsg ) then
!PUN          call msgerr(4, 'inconsistency found in SwanPunCollect: size of gathered data not correct ')
!PUN          return
!PUN       endif
!PUN    endif
!PUN    !
!PUN    ! calculate starting address of each local array with respect to the global array
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       idsplc(0) = 0
!PUN       do j = 1, MNPROC-1
!PUN          idsplc(j) = icount(j-1) + idsplc(j-1)
!PUN       enddo
!PUN    endif
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       allocate(    iarr(nvertsg))
!PUN       allocate(     arr(nvertsg))
!PUN       allocate(xcugrdgl(nvertsg))
!PUN       allocate(ycugrdgl(nvertsg))
!PUN    endif
!PUN    !
!PUN    ! gather different amounts of data from each processor to the master
!PUN    !
!PUN    call MPI_GATHERV( ivertp, nownv, MPI_INTEGER, iarr, icount, idsplc, MPI_INTEGER, 0, COMM, ierr )
!PUN    if ( ierr == MPI_SUCCESS ) call MPI_GATHERV( blknd, nownv, MPI_REAL, arr, icount, idsplc, MPI_REAL, 0, COMM, ierr )
!PUN    if ( ierr /= MPI_SUCCESS ) then
!PUN       write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
!PUN       call msgerr( 4, trim(msgstr) )
!PUN       return
!PUN    endif
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       do j = 1, nvertsg
!PUN          blkndc(iarr(j)) = arr(j)
!PUN       enddo
!PUN    endif
!PUN    !
!PUN    call MPI_GATHERV( xpl, nownv, MPI_REAL, arr, icount, idsplc, MPI_REAL, 0, COMM, ierr )
!PUN    if ( ierr /= MPI_SUCCESS ) then
!PUN       write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
!PUN       call msgerr( 4, trim(msgstr) )
!PUN       return
!PUN    endif
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       do j = 1, nvertsg
!PUN          xcugrdgl(iarr(j)) = arr(j)
!PUN       enddo
!PUN    endif
!PUN    !
!PUN    call MPI_GATHERV( ypl, nownv, MPI_REAL, arr, icount, idsplc, MPI_REAL, 0, COMM, ierr )
!PUN    if ( ierr /= MPI_SUCCESS ) then
!PUN       write (msgstr, '(a,i6)') ' MPI produces some internal error - return code is ',ierr
!PUN       call msgerr( 4, trim(msgstr) )
!PUN       return
!PUN    endif
!PUN    !
!PUN    if ( IAMMASTER ) then
!PUN       do j = 1, nvertsg
!PUN          ycugrdgl(iarr(j)) = arr(j)
!PUN       enddo
!PUN    endif
!PUN    !
!PUN    deallocate(ivertp,blknd)
!PUN    deallocate(xpl,ypl)
!PUN    if ( IAMMASTER ) deallocate(icount,idsplc,iarr,arr)
!PUN    !
end subroutine SwanPunCollect
