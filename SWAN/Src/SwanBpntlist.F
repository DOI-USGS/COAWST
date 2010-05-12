subroutine SwanBpntlist
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
!     Copyright (C) 2009  Delft University of Technology
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
!
!   Authors
!
!   40.80: Marcel Zijlema
!   40.92: Marcel Zijlema
!
!   Updates
!
!   40.80, April 2008: New subroutine
!   40.92,  June 2008: changes with respect to boundary polygons
!
!   Purpose
!
!   Makes list of boundary vertices in ascending order
!   - counterclockwise in case of sea/mainland boundaries
!   - clockwise in case of island boundaries
!
!   Method
!
!   The grid contains a number of boundary polygons
!   They are by definition closed
!   The first boundary polygon refers to sea/mainland boundary and the other polygons refers to island boundaries
!
!   The vertices which define the sea/mainland boundary are inserted in the counterclockwise direction
!   The vertices which define the island boundary are inserted in the clockwise direction
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Local variables
!
    integer                            :: icell      ! cell index
    integer, save                      :: ient = 0   ! number of entries in this subroutine
    integer                            :: iface      ! face index
    integer                            :: istat      ! indicate status of allocation
    integer                            :: j          ! loop counter
    integer                            :: k          ! counter
    integer, dimension(1)              :: kx         ! location of minimum value in array of x-coordinates of boundary vertices
    integer, dimension(1)              :: ky         ! location of minimum value in array of y-coordinates of boundary vertices
    integer                            :: m          ! loop counter
    integer                            :: maxnbp     ! maximum number of boundary vertices in set of polygons
    integer                            :: nbptot     ! total number of boundary vertices
    integer                            :: nptemp     ! auxiliary integer to store number of points temporarily
    integer, dimension(3)              :: v          ! vertices in present cell
    integer                            :: v1         ! first vertex of present face
    integer                            :: v2         ! second vertex of present face
    integer                            :: vc         ! considered vertex
    integer                            :: vcf        ! first considered vertex of a boundary polygon
    integer                            :: vn         ! next vertex with respect to considered vertex (counterclockwise)
    !
    integer, dimension(:), allocatable :: blistot    ! list of all boundary vertices in ascending order
    !
    real                               :: d1         ! distance of a point to origin
    real                               :: d2         ! distance of another point to origin
    !
    character(80)                      :: msgstr     ! string to pass message
    !
    logical                            :: firstvert  ! indicate whether considered vertex is first vertex of boundary polygon
    !
    type(celltype), dimension(:), pointer :: cell    ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face    ! datastructure for faces with their attributes
    type(verttype), dimension(:), pointer :: vert    ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanBpntlist')
    !
    ! if list of boundary vertices is already filled, return
    !
    if (allocated(blist)) return
    !
    ! point to vertex, cell and face objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    vert(:)%atti(BINDX) = 0
    vert(:)%atti(BPOL)  = 0
    nbpt                = 0
    !
    ! determine total number of boundary vertices
    !
    nbptot = count(mask=vert(:)%atti(VMARKER)==1)
    !
    allocate(blistot(nbptot))
    blistot = 0
    !
    ! determine first boundary vertex nearest to the origin
    !
    kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VMARKER)==1)
    ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VMARKER)==1)
    !
    if ( kx(1) == ky(1) ) then
       !
       vc = kx(1)
       !
    else
       !
       d1 = sqrt((vert(kx(1))%attr(VERTX))**2+(vert(kx(1))%attr(VERTY))**2)
       d2 = sqrt((vert(ky(1))%attr(VERTX))**2+(vert(ky(1))%attr(VERTY))**2)
       !
       if ( d1 < d2 ) then
          vc = kx(1)
       else
          vc = ky(1)
       endif
       !
    endif
    !
    ! store first boundary vertex
    !
    vcf                 = vc
    blistot(1)          = vc
    nbpol               = 1
    firstvert           = .true.
    vert(vc)%atti(BPOL) = nbpol
    !
    nptemp = 0
    !
    ! store next subsequent boundary vertices in ascending order
    !
    k     = 1
    iface = 1
    !
    faceloop: do
       !
       if ( face(iface)%atti(FMARKER) == 1 ) then
          !
          if ( firstvert ) then
             !
             icell = face(iface)%atti(FACEC1)
             !
             v(1) = cell(icell)%atti(CELLV1)
             v(2) = cell(icell)%atti(CELLV2)
             v(3) = cell(icell)%atti(CELLV3)
             !
             ! pick up next vertex (counterclockwise counting of vertices is assumed)
             !
             vn = 0
             do j = 1, cell(icell)%nov
                if ( v(j) == vc ) then
                   vn = v(mod(j,cell(icell)%nov)+1)
                   exit
                endif
             enddo
             !
             if ( vn == 0 ) goto 10
             if ( vert(vn)%atti(VMARKER) /= 1 ) goto 10
             firstvert = .false.
             !
          endif
          !
          v1 = face(iface)%atti(FACEV1)
          v2 = face(iface)%atti(FACEV2)
          !
          if ( v1 == vc ) then
             !
             if ( v2 == vcf ) vc = vcf
             if ( any( v2 == blistot ) ) goto 10
             !
             k = k + 1
             blistot(k)          = v2
             vert(v2)%atti(BPOL) = nbpol
             vc = v2
             !
          elseif ( v2 == vc ) then
             !
             if ( v1 == vcf ) vc = vcf
             if ( any( v1 == blistot ) ) goto 10
             !
             k = k + 1
             blistot(k)          = v1
             vert(v1)%atti(BPOL) = nbpol
             vc = v1
             !
          elseif ( vc == vcf ) then        ! end of considered boundary polygon is found
             !
             if ( any( v1 == blistot ) .and. any( v2 == blistot ) ) goto 10
             !
             ! store number of boundary vertices for present polygon
             !
             nbpt(nbpol) = k - nptemp
             nptemp = k
             !
             ! take first vertex of next boundary polygon
             !
             vc                  = v1
             vcf                 = vc
             k = k + 1
             blistot(k)          = vc
             nbpol               = nbpol + 1
             firstvert           = .true.
             vert(vc)%atti(BPOL) = nbpol
             !
             ! give error if more than 100 boundary polygons are found
             !
             if ( nbpol > 100 ) call msgerr ( 2, ' More than 100 boundary polygons are found in grid' )
             !
          endif
          !
          if ( k == nbptot ) exit faceloop
          !
       endif
       !
 10    continue
       iface = iface + 1
       if ( iface > nfaces ) iface = 1
       !
    enddo faceloop
    !
    ! store number of boundary vertices for last polygon
    !
    nbpt(nbpol) = nbptot - nptemp
    !
    ! check if list contains boundary vertices only
    !
    do j = 1, nbptot
       !
       vc = blistot(j)
       if (vert(vc)%atti(VMARKER) /= 1) then
          write (msgstr, '(a,i4,a)') ' Vertex with index ',vc,' in boundary list is not a valid boundary point'
          call msgerr( 2, trim(msgstr) )
       endif
       !
    enddo
    !
    ! determine maximum number of boundary vertices in set of polygons and allocate blist
    !
    maxnbp = maxval(nbpt)
    !
    if(.not.allocated(blist)) allocate (blist(maxnbp,nbpol), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanBpntlist: array blist ' )
       return
    endif
    blist = 0
    !
    ! fill blist in appropriate manner
    !
    k = 0
    !
    do j = 1, nbpol
       !
       do m = 1, nbpt(j)
          vc                   = blistot(k+m)
          blist(m,j)           = vc
          vert(vc)%atti(BINDX) = m
       enddo
       !
       k = k + nbpt(j)
       !
    enddo
    !
    deallocate(blistot)
    !
end subroutine SwanBpntlist
