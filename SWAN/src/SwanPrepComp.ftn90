subroutine SwanPrepComp ( cross )
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
!   Does some preparations before computation is started
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use m_obsta, only: OBSTDONE
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer, dimension(nfaces), intent(out) :: cross ! contains sequence number of obstacles for each face
                                                     ! where they crossing or zero if no crossing
!
!   Local variables
!
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: ivert    ! loop counter over vertices
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanPrepComp')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    ! deallocate arrays kvertc and kvertf (we don't use them anymore!)
    !
    if (allocated(kvertc)) deallocate(kvertc)
    if (allocated(kvertf)) deallocate(kvertf)
    !
    ! ghost and exception vertices are regarded as vertices with boundary condition
    !
    do ivert = 1, nverts
       if ( vmark(ivert) >= excmark ) vert(ivert)%atti(VBC) = 1
    enddo
    !
    ! find obstacles in computational grid, if present
    !
    if ( NUMOBS > 0 .and. .not.OBSTDONE ) then
       call SwanFindObstacles ( cross )
       OBSTDONE = .true.
    endif
    !
end subroutine SwanPrepComp
