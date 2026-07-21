subroutine SwanCheckGrid
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
!     Copyright (C) 1993-2024  Delft University of Technology
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program. If not, see <http://www.gnu.org/licenses/>.
!
!
!   Authors
!
!   40.80: Marcel Zijlema
!
!   Updates
!
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   Checks whether the mesh is suited for computation
!
!   Method
!
!   For the following aspects the triangular mesh is checked for:
!       1) the number of cells around vertex should be at least 4 and not larger than 10
!       2) the angles in each triangle should be smaller than 143 degrees
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Local variables
!
    integer                               :: i         ! loop counter
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: j         ! counter
    integer                               :: noc       ! number of cells around considered vertex
    integer                               :: v1        ! first vertex of present cell
    integer                               :: v2        ! second vertex of present cell
    integer                               :: v3        ! third vertex of present cell
    !
    real                                  :: cosphi1   ! cosine of the angle of first vertex in a triangle
    real                                  :: cosphi2   ! cosine of the angle of second vertex in a triangle
    real                                  :: cosphi3   ! cosine of the angle of third vertex in a triangle
    real                                  :: len12     ! squared length of face between vertices 1 and 2
    real                                  :: len13     ! squared length of face between vertices 1 and 3
    real                                  :: len23     ! squared length of face between vertices 2 and 3
    real                                  :: xdif12    ! difference in x-coordinate of vertices 1 and 2
    real                                  :: xdif13    ! difference in x-coordinate of vertices 1 and 3
    real                                  :: xdif23    ! difference in x-coordinate of vertices 2 and 3
    real                                  :: ydif12    ! difference in y-coordinate of vertices 1 and 2
    real                                  :: ydif13    ! difference in y-coordinate of vertices 1 and 3
    real                                  :: ydif23    ! difference in y-coordinate of vertices 2 and 3
    !
    logical                               :: acute     ! indicates whether triangles are acute (.TRUE.) or not (.FALSE.)
    logical                               :: badvertex ! indicates vertex has too less cells surrounded
    !
    character(80)                         :: msgstr    ! string to pass message
    !
    type(verttype), dimension(:), pointer :: vert   ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanCheckGrid')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    ! check whether the number of cells that meet at each internal vertex is at least 4
    ! (vertices at the boundaries not taken into account)
    ! check also whether the number of cells around each vertex is not larger than 10
    !
    badvertex = .false.
    do i = 1, nverts
       noc = vert(i)%noc
       if ( noc > 0 .and. noc < 4 .and. vmark(i) == 0 ) then
          if ( ITEST >= 30 .and. .not.badvertex ) write(PRINTF,'(a)')
          badvertex = .true.
          if ( ITEST >= 30 ) write(PRINTF,101) noc, i
       endif
       if ( noc > 10 ) then
          if ( ITEST >= 30 .and. .not.badvertex ) write(PRINTF,'(a)')
          badvertex = .true.
          if ( ITEST >= 30 ) write(PRINTF,101) noc, i
       endif
    enddo
    !
    if ( badvertex ) call msgerr (1, 'number of cells around vertex is smaller than 4 or larger than 10')
    !
    acute = .true.
    !
    ! check whether the angles in each triangle are smaller than 143 degrees
    !
    j = 0
    !
    do i = 1, ncells
       !
       v1 = kvertc(1,i)
       v2 = kvertc(2,i)
       v3 = kvertc(3,i)
       !
       xdif12 = xcugrd(v2) - xcugrd(v1)
       ydif12 = ycugrd(v2) - ycugrd(v1)
       xdif13 = xcugrd(v3) - xcugrd(v1)
       ydif13 = ycugrd(v3) - ycugrd(v1)
       xdif23 = xcugrd(v3) - xcugrd(v2)
       ydif23 = ycugrd(v3) - ycugrd(v2)
       !
       len12 = xdif12*xdif12 + ydif12*ydif12
       len13 = xdif13*xdif13 + ydif13*ydif13
       len23 = xdif23*xdif23 + ydif23*ydif23
       !
       ! is triangle acute ?
       !
       if (acute) acute = (len12+len23>len13) .and. (len23+len13>len12) .and. (len13+len12>len23)
       !
       cosphi1 =( xdif12*xdif13 + ydif12*ydif13)/(sqrt(len12*len13))
       cosphi2 =(-xdif12*xdif23 - ydif12*ydif23)/(sqrt(len12*len23))
       cosphi3 =( xdif13*xdif23 + ydif13*ydif23)/(sqrt(len13*len23))
       !
       if ( .not. cosphi1 > -0.8 .and. ITEST >= 100 ) then
          if ( j == 0 ) write(PRINTF,'(a)')
          write (msgstr,'(a,i6,a,f12.5,a)') 'angle of cell no. ',i,' is ',acos(cosphi1)/degrad,' degrees'
          call msgerr (1, trim(msgstr) )
          j = 1
       endif
       if ( .not. cosphi2 > -0.8 .and. ITEST >= 100 ) then
          if ( j == 0 ) write(PRINTF,'(a)')
          write (msgstr,'(a,i6,a,f12.5,a)') 'angle of cell no. ',i,' is ',acos(cosphi2)/degrad,' degrees'
          call msgerr (1, trim(msgstr) )
          j = 1
       endif
       if ( .not. cosphi3 > -0.8 .and. ITEST >= 100 ) then
          if ( j == 0 ) write(PRINTF,'(a)')
          write (msgstr,'(a,i6,a,f12.5,a)') 'angle of cell no. ',i,' is ',acos(cosphi3)/degrad,' degrees'
          call msgerr (1, trim(msgstr) )
          j = 1
       endif
       !
    enddo
    !
    if (acute) call msgerr (0, 'The grid contains solely acute triangles ')
    !
 101 format (' ',i2,' cells around vertex ',i6)
    !
end subroutine SwanCheckGrid
