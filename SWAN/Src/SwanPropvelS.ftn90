subroutine SwanPropvelS ( cad   , cas   , ux2   , uy2   , &
                          dep1  , dep2  , cax   , cay   , &
                          kwave , cgo   , spcsig, iddlow, &
                          iddtop, ecos  , esin  , coscos, &
                          sincos, sinsin, rdx   , rdy   , &
                          dhdx  , dhdy  , dkdx  , dkdy  )
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
!
!   Authors
!
!   40.59: Erick Rogers
!   40.80: Marcel Zijlema
!   41.02: Marcel Zijlema
!   41.06: Gerbrant van Vledder
!   41.07: Marcel Zijlema
!   41.35: Casey Dietrich
!   41.60: Marcel Zijlema
!
!   Updates
!
!   40.80,     July 2007: New subroutine
!   40.59,   August 2007: muddy bottom included
!   41.02, February 2009: adaption of velocities in case of diffraction
!   41.06,    March 2009: add option of limitation of velocity in theta-direction
!   41.07,   August 2009: add option of alternative formula for computation of
!                         wave transport velocity in theta-direction based on
!                         (x,y)-derivatives of the wave number
!                         (see Holthuijsen (2007), page 210, footnote 4)
!   41.35,    March 2012: add option of limitation on csigma and ctheta
!   41.60,     July 2015: more accurate computation of gradients of depth or wave number for turning rate
!
!   Purpose
!
!   computes wave transport velocities of energy in spectral space
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use m_diffr
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                        :: iddlow ! minimum direction bin that is propagated within a sweep
    integer, intent(in)                        :: iddtop ! maximum direction bin that is propagated within a sweep
    !
    real, dimension(MDC,MSC), intent(out)      :: cad    ! wave transport velocity in theta-direction
    real, dimension(MDC,MSC), intent(out)      :: cas    ! wave transport velocity in sigma-direction
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cax    ! wave transport velocity in x-direction
    real, dimension(MDC,MSC,ICMAX), intent(in) :: cay    ! wave transport velocity in y-direction
    real, dimension(MSC,ICMAX), intent(in)     :: cgo    ! group velocity
    real, dimension(MDC), intent(in)           :: coscos ! help array containing cosine to power 2 of spectral directions
    real, dimension(nverts), intent(in)        :: dep1   ! water depth at previous time level
    real, dimension(nverts), intent(in)        :: dep2   ! water depth at current time level
    real, intent(in)                           :: dhdx   ! derivative of depth in x-direction
    real, intent(in)                           :: dhdy   ! derivative of depth in y-direction
    real, dimension(MSC), intent(in)           :: dkdx   ! derivative of wave number in x-direction
    real, dimension(MSC), intent(in)           :: dkdy   ! derivative of wave number in y-direction
    real, dimension(MDC), intent(in)           :: ecos   ! help array containing cosine of spectral directions
    real, dimension(MDC), intent(in)           :: esin   ! help array containing sine of spectral directions
    real, dimension(MSC,ICMAX), intent(in)     :: kwave  ! wave number
    real, dimension(MDC), intent(in)           :: sincos ! help array containing sine * cosine of spectral directions
    real, dimension(MDC), intent(in)           :: sinsin ! help array containing sine to power 2 of spectral directions
    real, dimension(nverts), intent(in)        :: ux2    ! ambient velocity in x-direction at current time level
    real, dimension(nverts), intent(in)        :: uy2    ! ambient velocity in y-direction at current time level
    real, dimension(2), intent(in)             :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)             :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MSC), intent(in)           :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell    ! cell index
    integer                               :: id       ! loop counter over direction bins
    integer                               :: iddum    ! counter in directional space
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: iv1      ! first index in computational stencil
    integer                               :: iv2      ! second index in computational stencil
    integer                               :: iv3      ! third index in computational stencil
    !
    real                                  :: alpha    ! upper limit of CFL restriction
    real, dimension(3)                    :: cd       ! coefficients for computing cad
    real, dimension(10)                   :: cs       ! coefficients for computing cas
    real, dimension(3)                    :: dloc     ! local depth at vertices
    real                                  :: duxdx    ! derivative of ux2 to x
    real                                  :: duxdy    ! derivative of ux2 to y
    real                                  :: duydx    ! derivative of uy2 to x
    real                                  :: duydy    ! derivative of uy2 to y
    real                                  :: fac      ! a factor
    real                                  :: fac2     ! another factor
    real                                  :: frlim    ! frequency range in which limit on velocity in theta-direction is applied
    real                                  :: kd       ! help variable, wave number times water depth
    real                                  :: pp       ! power of the frequency dependent limiter on refraction
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanPropvelS')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! set velocities to zero
    !
    cad = 0.
    cas = 0.
    !
    ! if no frequency shift and no refraction, return
    !
    if ( (ITFRE == 0 .or. (.not.DYNDEP .and. ICUR == 0)) .and. IREFR == 0 ) return
    !
    ! initialize coefficients
    !
    cd = 0.
    cs = 0.
    !
    iv1 = vs(1)
    iv2 = vs(2)
    iv3 = vs(3)
    !
    dloc(1) = dep2(iv1)
    dloc(2) = dep2(iv2)
    dloc(3) = dep2(iv3)
    !
    ! if at least one vertex is dry, return
    !
    if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) return
    !
    ! compute the derivatives of the ambient current
    !
    if ( ICUR /= 0 ) then
       !
       duxdx = rdx(1) * (ux2(iv1) - ux2(iv2)) + rdx(2) * (ux2(iv1) - ux2(iv3))
       duxdy = rdy(1) * (ux2(iv1) - ux2(iv2)) + rdy(2) * (ux2(iv1) - ux2(iv3))
       duydx = rdx(1) * (uy2(iv1) - uy2(iv2)) + rdx(2) * (uy2(iv1) - uy2(iv3))
       duydy = rdy(1) * (uy2(iv1) - uy2(iv2)) + rdy(2) * (uy2(iv1) - uy2(iv3))
       !
    endif
    !
    ! compute wave transport velocity in sigma-direction
    !
    if ( ITFRE /= 0 .and. (DYNDEP .or. ICUR /= 0) ) then
       !
       ! compute time derivative of water depth
       !
       if ( DYNDEP ) cs(2) = ( dep2(iv1) - dep1(iv1) ) * RDTIM
       !
       ! compute coefficients depending on ambient currents
       !
       if ( ICUR /= 0 ) then
          !
          cs(3) = ux2(iv1) * dhdx
          cs(4) = uy2(iv1) * dhdy
          !
       endif
       !
       do is = 1, MSC
          !
          ! compute frequency-dependent coefficients
          !
          kd = min(30.,kwave(is,1) * dep2(iv1))
          !
          cs(1) =  kwave(is,1) * spcsig(is) / sinh (2.* kd)
          cs(5) = -cgo(is,1) * kwave(is,1)
          if ( IDIFFR /= 0 ) cs(5) = cs(5)*DIFPARAM(iv1)
          !
          cs( 6) = cs(1) * cs(2)
          cs( 7) = cs(1) * (cs(3)+cs(4))
          cs( 8) = cs(5) * duxdx
          cs( 9) = cs(5) * (duxdy+duydx)
          cs(10) = cs(5) * duydy
          !
          do iddum = iddlow-1, iddtop+1
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             cas(id,is) = cs(6)
             if ( ICUR /= 0 ) cas(id,is) = cas(id,is) + cs(7) + coscos(id)*cs(8) + sincos(id)*cs(9) + sinsin(id)*cs(10)
             !
          enddo
       enddo
       !
       ! limit velocity using Courant number
       !
       if ( int(PNUMS(33)) == 1 ) then
          !
          alpha = PNUMS(34)
          !
          do is = 1, MSC
             !
             fac2 = alpha * FRINTF * spcsig(is)
             !
             do iddum = iddlow-1, iddtop+1
                id = mod ( iddum - 1 + MDC , MDC ) + 1
                !
                fac = fac2 * ( abs((rdx(1)+rdx(2))*cax(id,is,1)) + abs((rdy(1)+rdy(2))*cay(id,is,1)) )
                !
                if ( abs(cas(id,is)) > fac ) cas(id,is) = cas(id,is) * fac / abs(cas(id,is))
                !
             enddo
             !
          enddo
          !
       endif
       !
    endif
    !
    ! compute wave transport velocity in theta-direction
    !
    if ( IREFR /= 0 ) then
       !
       do is = 1, MSC
          !
          ! compute frequency-dependent coefficients
          !
          if ( int(PNUMS(32)) == 0 ) then
             !
             kd = min(30.,kwave(is,1) * dep2(iv1))
             !
             cd(1) = spcsig(is) / sinh (2.* kd)
             !
             cd(2) = cd(1) * dhdx
             cd(3) = cd(1) * dhdy
             !
          else
             !
             cd(1) = -cgo(is,1) / kwave(is,1)
             !
             cd(2) = cd(1) * dkdx(is)
             cd(3) = cd(1) * dkdy(is)
             !
          endif
          !
          do iddum = iddlow-1, iddtop+1
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             cad(id,is) = esin(id)*cd(2) - ecos(id)*cd(3)
             if ( IDIFFR /= 0 ) cad(id,is) = cad(id,is)*DIFPARAM(iv1) - DIFPARDX(iv1)*cgo(is,1)*esin(id) + DIFPARDY(iv1)*cgo(is,1)*ecos(id)
             if ( ICUR   /= 0 ) cad(id,is) = cad(id,is) + sincos(id)*(duxdx-duydy) + sinsin(id)*duydx - coscos(id)*duxdy
             !
          enddo
          !
       enddo
       !
       ! adapt velocity in case of spherical coordinates
       !
       if ( KSPHER > 0 ) then
          !
          cd(1) = tan(DEGRAD*(vert(iv1)%attr(VERTY) + YOFFS)) / REARTH
          !
          do id = 1, MDC
             cd(2) = cd(1) * ecos(id)
             do is = 1, MSC
                cad(id,is) = cad(id,is) - cd(2)*(cax(id,is,1)*ecos(id) + cay(id,is,1)*esin(id))
             enddo
          enddo
          !
       endif
       !
       ! limit velocity in some frequency range if requested
       !
       if ( int(PNUMS(29)) == 1 ) then
          !
          frlim = PI2*PNUMS(26)
          pp    =     PNUMS(27)
          !
          do is = 1, MSC
             !
             fac = min(1.,(spcsig(is)/frlim)**pp)
             !
             do id = 1, MDC
                cad(id,is) = fac*cad(id,is)
             enddo
          enddo
          !
       endif
       !
       ! limit velocity using Courant number
       !
       if ( int(PNUMS(35)) == 1 ) then
          !
          alpha = PNUMS(36)
          !
          fac2 = alpha * DDIR
          !
          do is = 1, MSC
             !
             do iddum = iddlow-1, iddtop+1
                id = mod ( iddum - 1 + MDC , MDC ) + 1
                !
                fac = fac2 * ( abs((rdx(1)+rdx(2))*cax(id,is,1)) + abs((rdy(1)+rdy(2))*cay(id,is,1)) )
                !
                if ( abs(cad(id,is)) > fac ) cad(id,is) = cad(id,is) * fac / abs(cad(id,is))
                !
             enddo
             !
          enddo
          !
       endif
       !
    endif
    !
end subroutine SwanPropvelS
