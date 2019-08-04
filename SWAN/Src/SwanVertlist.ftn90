subroutine SwanVertlist ( compda )
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
!     Copyright (C) 1993-2019  Delft University of Technology
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
!   41.07: Casey Dietrich
!   41.48: Marcel Zijlema
!   41.68: Marcel Zijlema
!
!   Updates
!
!   40.80,   July 2007: New subroutine
!   41.07,   July 2009: small fix (assign ref.point to deepest point in case of no b.c.)
!   41.48,  March 2013: including order along a user-given direction
!   41.68, August 2015: introduction of a fixed number of sweeps per iteration
!   41.68,  April 2018: removal of the wavefront approach (based on reference point)
!
!   Purpose
!
!   Makes vertex list in line with sweep direction
!   Note: first sweep direction always equals user-given/wave/wind direction
!
!   Method
!
!   Sorting based on increasing distance along sweep direction
!
!   Modules used
!
    use ocpcomm4
    use swcomm2, only: COSWC, SINWC, VARWI
    use swcomm3, only: MCMVAR, JWX2, JWY2, JWX3, JWY3
    use m_genarr
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!PUN    use SIZES, only: MNPROC
!
    implicit none
!
!   Argument variables
!
    real, dimension(nverts,MCMVAR), intent(in) :: compda ! array containing space-dependent info (e.g. wind)
!
!   Local variables
!
    integer, save                     :: ient = 0   ! number of entries in this subroutine
    integer                           :: ierr       ! error indicator: ierr=0: no error, otherwise error
    integer                           :: istat      ! indicate status of allocation
    integer                           :: itmp       ! temporary stored integer for swapping
    integer                           :: j          ! loop counter over vertices
    integer                           :: k          ! counter
    integer, dimension(1)             :: kd         ! location of minimum value in array dist
    integer                           :: swpdir     ! sweep counter
    !
    real                              :: rtmp       ! temporary stored real for swapping
    real                              :: sdir       ! sweep direction
    real                              :: wdsum      ! total sum of wind direction
    real                              :: wx         ! wind velocity in x-direction
    real                              :: wy         ! wind velocity in y-direction
    !
    real, dimension(:,:), allocatable :: dist       ! distance of each point with respect to reference point
    !
    type(verttype), dimension(:), pointer :: vert   ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanVertlist')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    istat = 0
    if(.not.allocated(vlist)) allocate (vlist(nverts,nsweep), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanVertlist: array vlist ' )
       return
    endif
    !
    allocate (dist(nverts,nsweep))
    !
    ! check first sweep direction
    !
    if ( .not. asort > -999. ) then
!      if asort still does not have a value, try space-varying wind and take the mean of wind direction
       if ( VARWI ) then
!         retrieve current wind field
          if ( NSTATM == 1 ) call FLFILE ( 5, 6, WXI, WYI, 0, JWX2, JWX3, 0, JWY2, JWY3, COSWC, SINWC, compda, XCGRID, YCGRID, KGRPNT, ierr )
          k     = 0
          wdsum = 0.
          do j = 1, nverts
!            internal vertices and ghost vertices only
             if ( vmark(j) == 0 .or. vmark(j) == 999 ) then
                wx = compda(j,JWX2)
                wy = compda(j,JWY2)
                if ( wx /= 0. .or. wy /= 0. ) then
                   k = k + 1
                   wdsum = wdsum + atan2(wy,wx)
                endif
             endif
          enddo
          asort = wdsum / real(k)
!PUN          call SwanSumOverNodes ( asort )
!PUN          asort = asort / real(MNPROC)
       else
!         final attempt: set sweep direction to zero
          asort = 0.
       endif
    endif
    if ( ITEST >= 40 ) write (PRINTF,10) nsweep, 180.*asort/PI
    asort = asort - PI/real(nsweep)
    !
    ! order vertices according to sweep direction; base vector is user-given/wave/wind direction
    !
    sdir = asort + PI/real(nsweep)
    do swpdir = 1, nsweep
       do j = 1, nverts
          dist(j,swpdir) = vert(j)%attr(VERTX) * cos(sdir) + vert(j)%attr(VERTY) * sin(sdir)
       enddo
       sdir = sdir + PI2/real(nsweep)
    enddo
    !
    ! sort vertex list in order of increasing distance
    !
    do swpdir = 1, nsweep
       !
       do j = 1, nverts
          vlist(j,swpdir) = j
       enddo
       !
       do j = 1, nverts-1
          !
          kd = minloc(dist(j:nverts,swpdir))
          k  = kd(1) + j-1
          !
          if ( k /= j ) then
             !
             rtmp            = dist(j,swpdir)
             dist(j,swpdir)  = dist(k,swpdir)
             dist(k,swpdir)  = rtmp
             !
             itmp            = vlist(j,swpdir)
             vlist(j,swpdir) = vlist(k,swpdir)
             vlist(k,swpdir) = itmp
             !
          endif
          !
       enddo
       !
    enddo
    !
    deallocate(dist)
    !
 10 format (' Number of sweeps = ',i2,'; chosen wave direction for sweeping: ',f7.2,' degrees')
    !
end subroutine SwanVertlist
