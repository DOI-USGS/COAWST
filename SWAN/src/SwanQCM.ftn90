! This file contains data and routines for quasi-coherent modelling (QCM)
!
module SwanQCM
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!   41.90, June 2021: New module
!
!   Purpose
!
!   Contains data with respect to quasi-coherent modelling
!
!   Method
!
!   MODULE construct
!
!   Modules used
!
    implicit none
!
!   Module variables
!
    integer                                      :: mkxc      ! size of wave number space in x-direction
    integer                                      :: mkyc      ! size of wave number space in y-direction
    integer                                      :: mxd       ! size of dissipation grids in x-direction
    integer                                      :: myd       ! size of dissipation grids in y-direction
    integer                                      :: ncoz      ! size of coherent region / scattering wave number space
    !
    integer                                      :: lensav    ! size of work array wsave (FFT)
    integer                                      :: lensvd    ! size of work array wsavd (FFT)
    integer                                      :: lenwfd    ! size of work array wfd (FFT)
    integer                                      :: lenwft    ! size of work array wft (FFT)
    !
    real, dimension(:), save, allocatable        :: disbk0    ! bulk dissipation at previous iteration intended for surf breaking (global domain)
    real, dimension(:), save, allocatable        :: disbk1    ! bulk dissipation at current iteration intended for surf breaking (per subdomain)
    !
    real, dimension(:), save, allocatable        :: kx        ! x-component of wave number space
    real, dimension(:), save, allocatable        :: ky        ! y-component of wave number space
    !
    real, dimension(:), save, allocatable        :: xpsc      ! coherent scattering region for FFT (centered at active grid point)
    !
    real, dimension(:), save, allocatable        :: xpd       ! x-component of spatial lag grid intended for surf breaking
    real, dimension(:), save, allocatable        :: ypd       ! y-component of spatial lag grid intended for surf breaking
    !
    real, dimension(:), save, allocatable        :: kxd       ! x-component of wave number space asocciated with spatial lag grid
    real, dimension(:), save, allocatable        :: kyd       ! y-component of wave number space asocciated with spatial lag grid
!
!   Source text
!
contains
!
subroutine SWQCINIT ( BGRIDP, COMPDA )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!    41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!    41.91: Marcel Zijlema
!
!   Updates
!
!    41.90,     June 2021: New subroutine
!    41.91, February 2022: adding QC surf breaking
!
!   Purpose
!
!   Initializes quasi-coherent framework
!
!   Method
!
!   The coherent effects depend mainly on the nature of the medium variations, such as
!   depth and current (though still varying over many wave lengths), and the width of
!   the (incoming) spectrum
!
!   In general, narrow-banded wave groups in shallow water conditions cause O(1) variations
!   in the wave statistics over many wave lengths, implying that the wave field remains
!   correlated over such distances
!
!   To determine the importance of the coherent effects, knowledge of the spectral width
!   of the incident waves is required
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use m_genarr
    use m_parall
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, dimension(6*NBGRPT)    , intent(in) :: BGRIDP ! data concerning boundary grid points
    real   , dimension(MCGRD,MCMVAR), intent(in) :: COMPDA ! array containing space-dependent data
!
!   Local variables
!
    integer                         :: ib                  ! loop counter over boundary points
    integer                         :: ic                  ! computational grid point index
    integer                         :: id                  ! loop counter over directional bins
    integer, save                   :: ient = 0            ! number of entries in this subroutine
    integer                         :: ikx                 ! loop counter over x-components of wave number space
    integer                         :: iky                 ! loop counter over y-components of wave number space
    integer                         :: istat               ! indicate status of allocation
    integer                         :: is                  ! loop counter over frequency bins
    integer                         :: ix                  ! point index in x-direction
    integer                         :: ixd                 ! loop counter over x-components of spatial lag grid
    integer                         :: iy                  ! point index in y-direction
    integer                         :: iyd                 ! loop counter over y-components of spatial lag grid
    integer                         :: j                   ! loop counter
    integer                         :: n                   ! squared size coherent scattering region
    !
    real                            :: alfa                ! multiple of mean wave number to truncate scattering wave number space
    real, dimension(MSC)            :: arr                 ! auxiliary array
    real, dimension(MSC)            :: cg                  ! group velocity
    real                            :: clen                ! coherent length
    real                            :: dir1                ! first direction in directional sector
    real                            :: dir2                ! second direction in directional sector
    real                            :: dkx                 ! resolution of wave number space in x-direction
    real                            :: dky                 ! resolution of wave number space in y-direction
    real                            :: dp                  ! local depth
    real                            :: dxd                 ! resolution of dissipation grids in x-direction
    real                            :: dxps                ! spatial step of coherent scattering region (= pi/qmax)
    real                            :: dyd                 ! resolution of dissipation grids in y-direction
    real, dimension(MSC)            :: ebk                 ! variance density in wave number at boundary
    real                            :: fac                 ! auxiliary factor
    real                            :: hi                  ! incident wave height
    real                            :: hmax                ! maximum depth found
    real                            :: hmin                ! minimum depth found
    real                            :: hsi                 ! maximum incident wave height
    real, dimension(MSC)            :: k                   ! wave numbers
    real                            :: km                  ! local mean wave number
    real                            :: kmax                ! maximum wave number found
    real                            :: kmin                ! minimum wave number found
    real                            :: kmean               ! global mean wave number
    real                            :: kxlen               ! length of the wave number grid in x-direction
    real                            :: kxmax               ! upper bound of x-component in wave number space
    real                            :: kxmin               ! lower bound of x-component in wave number space
    real                            :: kylen               ! length of the wave number grid in y-direction
    real                            :: kymax               ! upper bound of y-component in wave number space
    real                            :: kymin               ! lower bound of y-component in wave number space
    real                            :: m0                  ! zeroth moment
    real                            :: qmax                ! integral bound of scattering wave number space
                                                           ! expressed as number of bins per direction
    real                            :: qmu                 ! maximum scattering wave number prescribed by user
    real                            :: rfac                ! resolution factor to determine step size of wave number grid
    real                            :: sd                  ! local standard deviation in wave number (1/m)
    real                            :: sdev                ! global standard deviation in wave number (1/m)
    !
    logical                         :: lpb                 ! indicate whether boundary point is a test point
    logical                         :: stpnow              ! indicate whether program must be terminated or not
    !
    character(80)                   :: msgstr              ! string to pass message
!
!   Structure
!
!   compute the mean wave number and standard deviation in wave number of the incoming spectrum
!   construct wave number space for scattering
!   compute coherent length scale
!   construct coherent scattering zone
!   construct grids for surf breaking
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCINIT')
    !
    ! compute the mean and standard deviation of the incoming spectrum in terms of wave number
    !
    kmean = -9999.
    sdev  =  9999.
    !
    ! also compute the maximum incident wave height meant for stopping criterion
    !
    hsi   = -9999.
    !
    if ( NBGRPT > 0 ) then
       !
       do ib = 1, NBGRPT
          !
          ic = BGRIDP(6*ib-5)
          !
          if ( BGRIDP(6*ib-4) == 1 ) then
             !
             ! first, compute moment of zeroth order of incoming spectrum
             !
             m0 = 0.
             do is = 1, MSC
                fac = SPCSIG(is)**2
                do id = 1, MDC
                   m0 = m0 + fac * AC2(id,is,ic)
                enddo
             enddo
             m0 = m0 * FRINTF * DDIR
             !
             if ( .not. m0 /= 0. ) cycle
             !
             ! compute the significant wave height
             !
             if ( m0 > 0. ) then
                hi = 4. * sqrt(m0)
             else
                hi = 0.
             endif
             !
             hsi = max ( hsi, hi )
             !
             ! next, obtain local depth ...
             !
             dp = COMPDA(ic,JDP2)
             !
             ! ... and subsequently the wave numbers
             !
             call KSCIP1 ( MSC, SPCSIG, dp, k, cg, arr, arr )
             !
             ! compute the incoming energy density in wave number space, integrated over all directions
             !
             do is = 1, MSC
                ebk(is) = 0.
                do id = 1, MDC-1
                   ebk(is) = ebk(is) + AC2(id,is,ic) + AC2(id+1,is,ic)
                enddo
                ! the Jacobian is the group velocity
                ebk(is) = 0.5 * DDIR * cg(is) * SPCSIG(is) * ebk(is)
             enddo
             !
             ! compute the mean wave number
             !
             km = 0.
             do is = 1, MSC-1
                km = km + ( k(is+1) - k(is) ) * ( k(is)*ebk(is) + k(is+1)*ebk(is+1) )
             enddo
             km = 0.5 * km / m0
             !
             kmean = max ( kmean, km )
             !
             ! compute the standard deviation
             !
             sd = 0.
             do is = 1, MSC-1
                sd = sd + ( k(is+1) - k(is) ) * ( (( k(is)- km )**2)*ebk(is) + (( k(is+1) - km )**2)*ebk(is+1) )
             enddo
             sd = sqrt( 0.5 * sd / m0 )
             !
             sdev = min ( sdev, sd )
             !
             ! test output: parameters in test points on boundary
             !
             if ( NPTST > 0 ) then
                do j = 1, NPTST
                   if ( OPTG /= 5 ) then
                      ix  = xytst(2*j-1)
                      iy  = xytst(2*j)
                      lpb = ic == KGRPNT(ix,iy)
                   else
                      ix  = xytst(j)
                      lpb = ic == ix
                   endif
                   if ( lpb ) then
                      if ( OPTG /= 5 ) then
                         write (PRTEST,101) ib, ix-1, iy-1
                      else
                         write (PRTEST,102) ib, ix
                      endif
                      write (PRTEST,103) km, sd
                   endif
                enddo
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! perform global reductions in parallel run
    call SWREDUCE( kmean, 1, SWREAL, SWMAX )
    call SWREDUCE( sdev , 1, SWREAL, SWMIN )
    call SWREDUCE( hsi  , 1, SWREAL, SWMAX )
    if ( stpnow() ) return
    !
    if ( sdev /= 9999. ) then
       if ( ITEST >= 50 ) write (PRTEST,'(a,f8.5)') ' standard deviation of incoming spectra (1/m) = ', sdev
    else
       call msgerr ( 3, 'unable to compute standard deviation of incoming spectra' )
       call msgerr ( 0, 'please instead specify resolution of the wave number grid!' )
    endif
    !
    ! redefine absolute stopping criterion by taking into account the scale of incident wave height
    !
    PNUMS(2) = PNUMS(2) * hsi
    !
    ! determine size and resolution of the wave number grid
    !
    kxmin = pscat(3)
    kymin = pscat(4)
    kxlen = pscat(5)
    kylen = pscat(6)
    !
    rfac  = pscat(7)
    !
    if ( .not. kxlen /= 0. .and. .not. kylen /= 0. ) then
       !
       ! construct wave number grid based on the frequency-direction grid
       !
       ! first, find minimum and maximum depth
       !
       hmin =  99999.
       hmax = -99999.
       !
       do ic = 1, MCGRD
          !
          dp = COMPDA(ic,JDP2)
          !
          if ( dp > DEPMIN ) then
             !
             hmin = min ( dp, hmin )
             hmax = max ( dp, hmax )
             !
          endif
          !
       enddo
       !
       ! perform global reductions in parallel run
       call SWREDUCE( hmin, 1, SWREAL, SWMIN )
       call SWREDUCE( hmax, 1, SWREAL, SWMAX )
       if ( stpnow() ) return
       !
       ! next, compute minimum and maximum absolute wave numbers ...
       !
       call KSCIP1 ( 1, spcsig(  1), hmin, k, cg, arr, arr )
       kmin = k(1)
       !
       call KSCIP1 ( 1, spcsig(MSC), hmax, k, cg, arr, arr )
       kmax = k(1)
       !
       ! finally, determine size and resolution of the intended wave number grid for scattering
       !
       if ( FULCIR ) then
          !
          kxmin = -kmax
          kxmax =  kmax
          kymin = -kmax
          kymax =  kmax
          !
       else
          !
          dir1 = spcdir(  1,1)
          dir2 = spcdir(MDC,1)
          !
          ! normalizes spectral directions to [0,2pi]
          !
          dir1 = dir1 - ( floor( dir1 / pi2 ) * pi2 )
          dir2 = dir2 - ( floor( dir2 / pi2 ) * pi2 )
          !
          if ( dir1 > 0. .and. .not. dir1 > 0.5*pi ) then
             !
             kxmax = kmax * cos( dir1 )
             kymax = kmax
             !
             if ( dir2 < pi ) then
                !
                kxmin = -kmax * cos( pi - dir2 )
                kymin = 0.
                !
             elseif ( dir2 < 1.5*pi ) then
                !
                kxmin = -kmax
                kymin = -kmax * cos( 1.5*pi - dir2 )
                !
             else
                !
                kxmax = max( kxmax, kmax * cos( dir2 ) )
                !
                kxmin = -kmax
                kymin = -kmax
                !
             endif
             !
          elseif ( dir1 > 0.5*pi .and. .not. dir1 > pi ) then
             !
             kxmin = -kmax
             kymax =  kmax * cos( 0.5*pi - dir1 )
             !
             if ( dir2 < 1.5*pi ) then
                !
                kxmax = 0.
                kymin = -kmax * cos( 1.5*pi - dir2 )
                !
             elseif ( dir2 < pi2 ) then
                !
                kxmax =  kmax * cos( dir2 )
                kymin = -kmax
                !
             else
                !
                kymax = max( kymax, kmax * sin( dir2 ) )
                !
                kxmax =  kmax
                kymin = -kmax
                !
             endif
             !
          elseif ( dir1 > pi .and. .not. dir1 > 1.5*pi ) then
             !
             kxmin = -kmax * cos( pi - dir1 )
             kymin = -kmax
             !
             if ( dir2 < pi2 ) then
                !
                kxmax = kmax * cos( dir2 )
                kymax = 0.
                !
             elseif ( dir2 < 0.5*pi ) then
                !
                kxmax = kmax
                kymax = kmax * sin( dir2 )
                !
             else
                !
                kxmin = min( kxmin, -kmax * cos( pi - dir2 ) )
                !
                kxmax = kmax
                kymax = kmax
                !
             endif
             !
          elseif ( dir1 > 1.5*pi .and. .not. dir1 > pi2 ) then
             !
             kxmax =  kmax
             kymin = -kmax * cos( 1.5*pi - dir1 )
             !
             if ( dir2 < 0.5*pi ) then
                !
                kxmin = 0.
                kymax = kmax * sin( dir2 )
                !
             elseif ( dir2 < pi ) then
                !
                kxmin = -kmax * cos( pi - dir2 )
                kymax =  kmax
                !
             else
                !
                kymin = min( kymin, -kmax * cos( 1.5*pi - dir2 ) )
                !
                kxmin = -kmax
                kymax =  kmax
                !
             endif
             !
          endif
          !
       endif
       !
       kxlen = kxmax - kxmin
       kylen = kymax - kymin
       !
    endif
    !
    if ( mkxc == -1 .and. mkyc == -1 ) then
       dkx  = sdev / rfac
       dky  = sdev / rfac
       mkxc = 1 + nint(kxlen/dkx)
       mkyc = 1 + nint(kylen/dky)
    else
       dkx  = kxlen / float(mkxc)
       dky  = kylen / float(mkyc)
       mkxc = mkxc + 1
       mkyc = mkyc + 1
    endif
    if ( dkx > sdev .or. dky > sdev ) then
       call msgerr (1,'spectral space too coarse to capture wave interferences properly')
       write (PRINTF, 104) PI2 / max(dkx,dky)
    endif
    !
    ! allocate wave number grid for scattering
    !
    istat = 0
    if (                  .not.allocated(kx) ) allocate(kx(mkxc), stat = istat)
    if ( istat == 0 .and. .not.allocated(ky) ) allocate(ky(mkyc), stat = istat)
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: wave number space and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
    ! setup the wave number grid
    !
    kx(1) = kxmin
    do ikx = 2, mkxc
       kx(ikx) = kx(ikx-1) + dkx
    enddo
    !
    ky(1) = kymin
    do iky = 2, mkyc
       ky(iky) = ky(iky-1) + dky
    enddo
    !
    ! compute coherent radius
    !
    clen = PI2 / max(dkx,dky)
    !
    if ( ITEST >= 50 ) write (PRTEST,'(a,f9.2)') ' maximum resolvable coherent length scale (m) = ',clen
    !
    ! construct coherent scattering region
    ! note: medium variations only affect statistics within a radius clen/2 around grid point (see Smit et al, 2015)
    !
    alfa = pscat(1)
    qmu  = pscat(2)
    !
    ! QC domain delimited by a multiple of the mean wave number or user-prescribed maximum scattering wave number
    qmax = min (qmu, alfa * kmean) / max(dkx,dky)
    if ( qmax < 2. ) then
       call msgerr (1,'range of the QC convolution integral too small - set to minimum required range')
       qmax = 2.
    endif
    ncoz = 2*floor(ceiling(qmax)/2.)   ! must be even number to execute the QC summation correctly (see routine SWQCSCAT)
    dxps = clen / real(ncoz) / 2.
    !
    if ( .not.allocated(xpsc) ) allocate( xpsc(ncoz) )
    !
    ! coherent zone is shifted with the origin centered
    !
    xpsc(1) = -clen / 4.
    do j = 2, ncoz
       xpsc(j) = xpsc(j-1) + dxps
    enddo
    !
    if (SCREEN /= PRINTF ) then
       write (PRINTF,105) clen, ncoz, xpsc(1), xpsc(ncoz)+dxps, xpsc(1), xpsc(ncoz)+dxps, -PI/dxps, PI/dxps, -PI/dxps, PI/dxps
       if (IAMMASTER .and. IGEN == 4 ) write (SCREEN,106) clen, ncoz
    endif
    !
    ! set up grids for surf breaking
    !
    if ( ISURF > 0 .and. IGEN == 4 ) then
       !
       mxd = 2*floor( (kxmin+kxlen)/dkx )
       myd = 2*floor( (kymin+kylen)/dky )
       !
       dxd = pi2 / ( mxd * dkx )
       dyd = pi2 / ( myd * dky )
       !
       ! allocate grids
       !
       istat = 0
       if (                  .not.allocated(xpd) ) allocate(xpd(mxd+1), stat = istat)
       if ( istat == 0 .and. .not.allocated(ypd) ) allocate(ypd(myd+1), stat = istat)
       if ( istat == 0 .and. .not.allocated(kxd) ) allocate(kxd(mxd+1), stat = istat)
       if ( istat == 0 .and. .not.allocated(kyd) ) allocate(kyd(myd+1), stat = istat)
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: surf breaking grids and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       ! construct grids
       !
       do ixd = 1, mxd+1
          !
          xpd(ixd) = ( ixd-1 - mxd/2 ) * dxd
          kxd(ixd) = ( ixd-1 - mxd/2 ) * dkx
          !
       enddo
       !
       do iyd = 1, myd+1
          !
          ypd(iyd) = ( iyd-1 - myd/2 ) * dyd
          kyd(iyd) = ( iyd-1 - myd/2 ) * dky
          !
       enddo
       !
    else
       !
       mxd = 0
       myd = 0
       !
    endif
    !
    mxd = mxd + 1
    myd = myd + 1
    !
    ! allocate and initialize arrays for storing bulk dissipation
    !
    if ( ISURF > 0 .and. IGEN == 4 ) then
       !
       istat = 0
       if (                  .not.allocated(disbk0) ) allocate(disbk0(MCGRDGL), stat = istat)
       if ( istat == 0 .and. .not.allocated(disbk1) ) allocate(disbk1(MCGRD  ), stat = istat)
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: bulk dissipation and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       disbk0 = 0.
       disbk1 = 0.
       !
    endif
    !
    ! set size of work arrays for FFT related to scattering ...
    !
    n = ncoz * ncoz
    !
    lenwft = 2 * n
    lensav = 2 * ncoz + int( log(float(ncoz)) / log(2.) ) + 4
    lensav = 2 * lensav
    !
    ! ... and for FFT related to surf breaking
    !
    lenwfd = 2 * mxd * myd
    lensvd = 2 * mxd + int( log(float(mxd)) / log(2.) ) + 2 * myd + int( log(float(myd)) / log(2.) ) + 8
    !
    ! format statements
    !
 101 format (' boundary point', 3i8)
 102 format (' boundary vertex', 2i8)
 103 format (' mean wave number and standard deviation: ', 2e12.4)
 104 format (' maximum resolvable coherent length scale is ', f9.2,' m')
 105 format (/' preparing data for quasi-coherent modelling:'/           &
             ' - maximum resolvable coherent length scale = ',f9.2,' m'/ &
             ' - number of Fourier components             = ', i4/       &
             ' -  xmin = ',f9.2,                                         &
             ',  xmax = ',f9.2,                                          &
             ',  ymin = ',f9.2,                                          &
             ',  ymax = ',f9.2,' (m)'/                                   &
             ' - qxmin = ',f9.3,                                         &
             ', qxmax = ',f9.3,                                          &
             ', qymin = ',f9.3,                                          &
             ', qymax = ',f9.3,' (rad/m)')
 106 format (' ... preparing data for quasi-coherent modelling'//          &
             '     maximum resolvable coherent length scale = ',f9.2,' m'/ &
             '     number of Fourier components             = ', i4/)
    !
end subroutine SWQCINIT
!
subroutine SWQCDFT ( sigft, cgft, dep2, kwave, cgo, cft, rft, sft, wft, wsave )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!    41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!    41.90, June 2021: New subroutine
!
!   Purpose
!
!   Performs discrete Fourier transforms of modulation of intrinsic frequency and group velocity
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use m_genarr
    use SwanGriddata
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    complex(kind=8), dimension(ncoz,ncoz)    , intent(inout) :: cft   ! Fourier coefficients (FFT)
    real           , dimension(ncoz,ncoz,MSC), intent(out)   :: cgft  ! Fourier-transformed modulation of group velocity
    real           , dimension(MSC)          , intent(in)    :: cgo   ! group velocity in frequency-direction space
    real           , dimension(MCGRD)        , intent(in)    :: dep2  ! water depth at current time level
    real           , dimension(MSC)          , intent(in)    :: kwave ! wave number in frequency-direction space
    real(kind=8)   , dimension(ncoz,ncoz)    , intent(inout) :: rft   ! input data (FFT)
    real(kind=8)   , dimension(ncoz,ncoz)    , intent(inout) :: sft   ! input data (FFT)
    real           , dimension(ncoz,ncoz,MSC), intent(out)   :: sigft ! Fourier-transformed modulation of intrinsic frequency
    real(kind=8)   , dimension(lenwft)       , intent(inout) :: wft   ! work array (FFT)
    real(kind=8)   , dimension(lensav)       , intent(inout) :: wsave ! work array (FFT)
!
!   Local variables
!
    integer                        :: i        ! loop counter
    integer, save                  :: ient = 0 ! number of entries in this subroutine
    integer                        :: ierr     ! error indicator
    integer                        :: ikx      ! loop counter over x-components of wave number space
    integer                        :: iky      ! loop counter over y-components of wave number space
    integer                        :: is       ! loop counter over frequency bins
    integer                        :: ish      ! shift level for shifting matrix
                                               ! note: ish+1 is the centre of coherent scattering zone
    integer                        :: ix       ! point index in x-direction
    integer                        :: iy       ! point index in y-direction
    integer                        :: j        ! loop counter
    !
    real                           :: cg       ! group velocity within coherent zone
    real                           :: cgp      ! group velocity at current grid point
    real                           :: d        ! depth value  within coherent zone obtained from input grid
    real                           :: dp       ! depth at current grid point
    real                           :: hmax     ! maximum depth in coherent zone
    real                           :: hmin     ! minimum depth in coherent zone
    real                           :: k        ! wave number
    real                           :: sig      ! relative frequency within coherent zone
    real                           :: sigp     ! relative frequency at current grid point
    real                           :: SVALQI   ! function giving interpolated value of an input array
    real                           :: x        ! x-coordinate of point in coherent region
    real                           :: xp       ! x-coordinate of computational grid point
    real                           :: y        ! y-coordinate of point in coherent region
    real                           :: yp       ! y-coordinate of computational grid point
    real                           :: wl       ! water level within coherent zone obtained from input grid
    !
    character(80)                  :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCDFT')
    !
    sigft = 0.
    cgft  = 0.
    !
    if ( optg /= 5 ) then
       ix = IXCGRD(1)
       iy = IYCGRD(1)
       xp = xcgrid(ix,iy)
       yp = ycgrid(ix,iy)
    else
       ix = 0
       iy = 0
       xp = xcugrd(vs(1))
       yp = ycugrd(vs(1))
    endif
    !
    ierr = 0
    ish  = ncoz/2
    !
    dp = dep2(KCGRD(1))
    !
    if ( dp > DEPMIN ) then
       !
       do is = 1, MSC
          !
          k = kwave(is)
          !
          sigp = spcsig(is)
          cgp  = cgo(is)
          !
          hmin =  99999.
          hmax = -99999.
          !
          ! compute modulation of relative frequency and group velocity in coherent zone
          !
          rft = 0.
          sft = 0.
          !
          do i = 1, ncoz
             do j = 1, ncoz
                !
                x = xp + xpsc(i)
                y = yp + xpsc(j)
                !
                d = SVALQI ( x, y, 1, DEPTH, 1, ix, iy )
                !
                if ( d /= EXCFLD(1) ) then
                   !
                   if ( VARWLV ) then
                      wl = SVALQI ( x, y, 7, WLEVL, 1, ix, iy )
                      d  = d + wl
                   endif
                   ! add constant water level
                   d = d + WLEV
                   !
                   hmin = min ( d, hmin )
                   hmax = max ( d, hmax )
                   !
                   sig = sqrt( k * GRAV * tanh(k*d) )
                   cg  = ( sig / k ) * 0.5 * ( 1. + 2.*k*d / ( sinh(2.*k*d) ) )
                   !
                   ! note: sigp = sig(ish+1,ish+1) = spcsig; likewise, cgp = cg(ish+1,ish+1) = cgo
                   rft(i,j) = sig - sigp
                   sft(i,j) = cg  - cgp
                   !
                endif
                !
             enddo
          enddo
          !
          if ( k*hmin > 2.3 ) cycle
          !
          if ( (hmax-hmin)/(xpsc(ncoz)-xpsc(1)) < 0.001 ) cycle
          !
          ! to avoid discontinuities at boundaries of coherent zone apply Tukey window to relative frequency
          !
          call tukeywin ( rft, ncoz )
          !
          ! compute Fourier components
          !
          cft = cmplx(rft)
          !
          call cfft2f ( ncoz, ncoz, ncoz, cft, wsave, lensav, wft, lenwft, ierr )
          if ( ierr /= 0 ) goto 10
          !
          ! compute spectrum
          !
          cft = cshift(cft, shift=ish, dim=2 )
          cft = cshift(cft, shift=ish, dim=1 )
          !
          ! the phase of spectrum is altered owing to the shift of coherent zone
          !
          do i = 1, ncoz
             do j = 1, ncoz
                cft(i,j) = (-1.)**(i+j) * cft(i,j)
             enddo
          enddo
          !
          ! save only the imaginary part (to be used in the QC scatterer term)
          !
          sigft(:,:,is) = imag(cft(:,:))
          !
          ! to avoid discontinuities at boundaries of coherent zone apply Tukey window to group velocity
          !
          call tukeywin ( sft, ncoz )
          !
          ! compute Fourier components
          !
          cft = cmplx(sft)
          !
          call cfft2f ( ncoz, ncoz, ncoz, cft, wsave, lensav, wft, lenwft, ierr )
          if ( ierr /= 0 ) goto 10
          !
          ! compute spectrum
          !
          cft = cshift(cft, shift=ish, dim=2 )
          cft = cshift(cft, shift=ish, dim=1 )
          !
          ! the phase of spectrum is altered owing to the shift of coherent zone
          !
          do i = 1, ncoz
             do j = 1, ncoz
                cft(i,j) = (-1.)**(i+j) * cft(i,j)
             enddo
          enddo
          !
          ! save only the real part (to be used in the QC scatterer term)
          !
          cgft(:,:,is) = real(cft(:,:))
          !
       enddo
       !
    endif
    !
 10 if ( ierr /= 0 ) then
       write (msgstr, '(a,i6)') 'something went wrong with the Fourier transform - return code is ',ierr
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
end subroutine SWQCDFT
!
subroutine SWQCUFT ( uxft, uyft, dep2, ux2, uy2, cft, rft, sft, wft, wsave )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!    41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!    41.90, June 2021: New subroutine
!
!   Purpose
!
!   Performs discrete Fourier transforms of modulation of ambient currents
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use m_genarr
    use SwanGriddata
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    complex(kind=8), dimension(ncoz,ncoz), intent(inout) :: cft   ! Fourier coefficients (FFT)
    real           , dimension(MCGRD)    , intent(in)    :: dep2  ! water depth at current time level
    real   (kind=8), dimension(ncoz,ncoz), intent(inout) :: rft   ! input data (FFT)
    real   (kind=8), dimension(ncoz,ncoz), intent(inout) :: sft   ! input data (FFT)
    real           , dimension(MCGRD)    , intent(in)    :: ux2   ! ambient velocity in x-direction at current time level
    complex        , dimension(ncoz,ncoz), intent(out)   :: uxft  ! u-component of Fourier-transformed modulation of ambient current
    real           , dimension(MCGRD)    , intent(in)    :: uy2   ! ambient velocity in y-direction at current time level
    complex        , dimension(ncoz,ncoz), intent(out)   :: uyft  ! v-component of Fourier-transformed modulation of ambient current
    real   (kind=8), dimension(lenwft)   , intent(inout) :: wft   ! work array (FFT)
    real   (kind=8), dimension(lensav)   , intent(inout) :: wsave ! work array (FFT)
!
!   Local variables
!
    integer                        :: i        ! loop counter
    integer, save                  :: ient = 0 ! number of entries in this subroutine
    integer                        :: ierr     ! error indicator
    integer                        :: ish      ! shift level for shifting matrix
    integer                        :: ix       ! point index in x-direction
    integer                        :: iy       ! point index in y-direction
    integer                        :: j        ! loop counter
    !
    real                           :: cgmx     ! max Froude number times wave celerity
    real                           :: dp       ! local depth
    real                           :: fac      ! auxiliary factor
    real                           :: SVALQI   ! function giving interpolated value of an input array
    real                           :: u        ! u-component of ambient current within coherent zone obtained from input grid
    real                           :: uu       ! u-component of ambient current w.r.t. user coordinates
    real                           :: uxp      ! u-component of ambient current at current grid point
    real                           :: uyp      ! v-component of ambient current at current grid point
    real                           :: v        ! v-component of ambient current within coherent zone obtained from input grid
    real                           :: vtot     ! velocity magnitude
    real                           :: vv       ! v-component of ambient current w.r.t. user coordinates
    real                           :: x        ! x-coordinate of point in coherent region
    real                           :: xp       ! x-coordinate of computational grid point
    real                           :: y        ! y-coordinate of point in coherent region
    real                           :: yp       ! y-coordinate of computational grid point
    !
    character(80)                  :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCUFT')
    !
    uxft = (0.,0.)
    uyft = (0.,0.)
    !
    if ( optg /= 5 ) then
       ix = IXCGRD(1)
       iy = IYCGRD(1)
       xp = xcgrid(ix,iy)
       yp = ycgrid(ix,iy)
    else
       ix = 0
       iy = 0
       xp = xcugrd(vs(1))
       yp = ycugrd(vs(1))
    endif
    !
    ierr = 0
    ish  = ncoz/2
    !
    dp  = dep2(KCGRD(1))
    uxp = ux2 (KCGRD(1))
    uyp = uy2 (KCGRD(1))
    !
    cgmx = PNUMS(18) * sqrt( GRAV*dp )
    !
    if ( dp > DEPMIN ) then
       !
       ! compute modulation of currents in coherent zone
       !
       rft = 0.
       sft = 0.
       !
       do i = 1, ncoz
          do j = 1, ncoz
             !
             x = xp + xpsc(i)
             y = yp + xpsc(j)
             !
             u = SVALQI ( x, y, 2, UXB, 0, ix, iy )
             v = SVALQI ( x, y, 3, UYB, 0, ix, iy )
             !
             if ( .not. u /= 0. .and. .not. v /= 0. ) cycle
             !
             uu =  u*COSVC + v*SINVC
             vv = -u*SINVC + v*COSVC
             !
             vtot = sqrt( uu*uu + vv*vv )
             !
             if ( vtot > cgmx ) then
                fac = cgmx / vtot
                uu = fac * uu
                vv = fac * vv
             endif
             !
             rft(i,j) = uu - uxp
             sft(i,j) = vv - uyp
             !
          enddo
       enddo
       !
       ! to avoid discontinuities at boundaries of coherent zone apply Tukey window to u-component
       !
       call tukeywin ( rft, ncoz )
       !
       ! compute Fourier components
       !
       cft = cmplx(rft)
       !
       call cfft2f ( ncoz, ncoz, ncoz, cft, wsave, lensav, wft, lenwft, ierr )
       if ( ierr /= 0 ) goto 10
       !
       ! compute spectrum
       !
       cft = cshift(cft, shift=ish, dim=2 )
       cft = cshift(cft, shift=ish, dim=1 )
       !
       ! the phase of spectrum is altered owing to the shift of coherent zone
       !
       do i = 1, ncoz
          do j = 1, ncoz
             cft(i,j) = (-1.)**(i+j) * cft(i,j)
          enddo
       enddo
       !
       uxft = cft
       !
       ! to avoid discontinuities at boundaries of coherent zone apply Tukey window to v-component
       !
       call tukeywin ( sft, ncoz )
       !
       ! compute Fourier components
       !
       cft = cmplx(sft)
       !
       call cfft2f ( ncoz, ncoz, ncoz, cft, wsave, lensav, wft, lenwft, ierr )
       if ( ierr /= 0 ) goto 10
       !
       ! compute spectrum
       !
       cft = cshift(cft, shift=ish, dim=2 )
       cft = cshift(cft, shift=ish, dim=1 )
       !
       ! the phase of spectrum is altered owing to the shift of coherent zone
       !
       do i = 1, ncoz
          do j = 1, ncoz
             cft(i,j) = (-1.)**(i+j) * cft(i,j)
          enddo
       enddo
       !
       uyft = cft
       !
    endif
    !
 10 if ( ierr /= 0 ) then
       write (msgstr, '(a,i6)') 'something went wrong with the Fourier transform - return code is ',ierr
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    !
end subroutine SWQCUFT
!
subroutine QCSOURCE ( imatra, imatda, iter  , ac2   , dep2  , ux2   , uy2   , &
                      swpdir, ix    , iy    , rdx   , rdy   , kwave , cgo   , &
                      sigft , cgft  , uxft  , uyft  , memqcm, memqcb, plqcs , &
                      plwbrk, dissc0, dissc1, genc0 , genc1 , redc0 , redc1 , &
                      spcsig, spcdir, idcmin, idcmax, isstop, ecos  , esin  , &
                      etot  , hm    , qb    , smebrk, kteta , kmespc, cft   , &
                      rft   , sft   , wft   , wsave , cfd   , wfd   , wsavd   &
                                                                            )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!   41.91: Marcel Zijlema
!
!   Updates
!
!   41.90, October  2021: New subroutine
!   41.91, February 2022: extension QCM with surf breaking
!
!   Purpose
!
!   Computes the source terms in the quasi-coherent framework
!
!   Method
!
!   compute interactions due to medium variations and surf breaking
!
!   in the near future, other source terms, e.g. bed friction
!   and a nonlinear forcing term (triad), will be included
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
!
    implicit none
!
!   Argument variables
!
    integer                                  , intent(in   ) :: isstop ! maximum frequency that is propagated within a sweep
    integer                                  , intent(in   ) :: iter   ! iteration counter
    integer                                  , intent(in   ) :: ix     ! counter of grid points in x-direction
    integer                                  , intent(in   ) :: iy     ! counter of grid points in y-direction
    integer                                  , intent(in   ) :: swpdir ! sweep counter
    !
    integer        , dimension(MSC)          , intent(in   ) :: idcmax ! maximum frequency-dependent counter in directional space
    integer        , dimension(MSC)          , intent(in   ) :: idcmin ! minimum frequency-dependent counter in directional space
    !
    real                                     , intent(in   ) :: etot   ! total wave energy density
    real                                     , intent(in   ) :: hm     ! maximum wave height
    real                                     , intent(in   ) :: kmespc ! mean average wavenumber according to the WAM formulation
    real                                     , intent(in   ) :: kteta  ! number of directional partitions
    real                                     , intent(in   ) :: qb     ! fraction of breaking waves
    real                                     , intent(in   ) :: smebrk ! mean frequency according to first order moment
    !
    real           , dimension(MDC,MSC,MCGRD), intent(in   ) :: ac2    ! action density at current time level
    complex(kind=8), dimension(myd,mxd)      , intent(inout) :: cfd    ! Fourier coefficients (FFT) for surf breaking
    complex(kind=8), dimension(ncoz,ncoz)    , intent(inout) :: cft    ! Fourier coefficients (FFT) for scattering
    real           , dimension(ncoz,ncoz,MSC), intent(inout) :: cgft   ! Fourier-transformed modulation of group velocity
    real           , dimension(MSC,MICMAX)   , intent(in   ) :: cgo    ! group velocity
    real           , dimension(MCGRD)        , intent(in)    :: dep2   ! water depth at current time level
    real           , dimension(MDC,MSC,MDISP), intent(out  ) :: dissc0 ! explicit part of dissipation in present geographical point for output purposes
    real           , dimension(MDC,MSC,MDISP), intent(out  ) :: dissc1 ! implicit part of dissipation in present geographical point for output purposes
    real           , dimension(MDC)          , intent(in   ) :: ecos   ! help array containing cosine of spectral directions
    real           , dimension(MDC)          , intent(in   ) :: esin   ! help array containing sine of spectral directions
    real           , dimension(MDC,MSC,MGENR), intent(out  ) :: genc0  ! explicit part of generation in present geographical point for output purposes
    real           , dimension(MDC,MSC,MGENR), intent(out  ) :: genc1  ! implicit part of generation in present geographical point for output purposes
    real           , dimension(MSC,MICMAX)   , intent(in   ) :: kwave  ! wave number
    real           , dimension(MDC,MSC)      , intent(out  ) :: imatra ! coefficients of right hand side
    real           , dimension(MDC,MSC)      , intent(out  ) :: imatda ! coefficients of main diagonal of matrix
    real           , dimension(MDC,MSC,MCGRD), intent(inout) :: memqcb ! auxiliary array to store results of QC surf breaking in full spectral space
    real           , dimension(MDC,MSC,MCGRD), intent(inout) :: memqcm ! auxiliary array to store results of QC scattering in full spectral space
    real           , dimension(MDC,MSC,NPTST), intent(out  ) :: plqcs  ! array containing the QC scattering source term for test output
    real           , dimension(MDC,MSC,NPTST), intent(out  ) :: plwbrk ! array containing the surf breaking source term for test output
    real           , dimension(2)            , intent(in   ) :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real           , dimension(2)            , intent(in   ) :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real           , dimension(MDC,MSC,MREDS), intent(out  ) :: redc0  ! explicit part of redistribution in present geographical point for output purposes
    real           , dimension(MDC,MSC,MREDS), intent(out  ) :: redc1  ! implicit part of redistribution in present geographical point for output purposes
    real   (kind=8), dimension(ncoz,ncoz)    , intent(inout) :: rft    ! input data (FFT)
    real           , dimension(ncoz,ncoz,MSC), intent(inout) :: sigft  ! Fourier-transformed modulation of intrinsic frequency
    real           , dimension(MDC,6)        , intent(in   ) :: spcdir ! (*,1): spectral direction bins (radians)
                                                                       ! (*,2): cosine of spectral directions
                                                                       ! (*,3): sine of spectral directions
                                                                       ! (*,4): cosine^2 of spectral directions
                                                                       ! (*,5): cosine*sine of spectral directions
                                                                       ! (*,6): sine^2 of spectral directions
    real   (kind=8), dimension(ncoz,ncoz)    , intent(inout) :: sft    ! input data (FFT)
    real           , dimension(MSC)          , intent(in   ) :: spcsig ! relative frequency bins
    real           , dimension(MCGRD)        , intent(in   ) :: ux2    ! ambient velocity in x-direction at current time level
    complex        , dimension(ncoz,ncoz)    , intent(inout) :: uxft   ! u-component of Fourier-transformed modulation of ambient current
    real           , dimension(MCGRD)        , intent(in   ) :: uy2    ! ambient velocity in y-direction at current time level
    complex        , dimension(ncoz,ncoz)    , intent(inout) :: uyft   ! v-component of Fourier-transformed modulation of ambient current
    real   (kind=8), dimension(lenwfd)       , intent(inout) :: wfd    ! work array (FFT) for surf breaking
    real   (kind=8), dimension(lenwft)       , intent(inout) :: wft    ! work array (FFT) for scattering
    real   (kind=8), dimension(lensvd)       , intent(inout) :: wsavd  ! work array (FFT) for surf breaking
    real   (kind=8), dimension(lensav)       , intent(inout) :: wsave  ! work array (FFT) for scattering
!
!   Local variables
!
    integer, save                :: ient = 0 ! number of entries in this subroutine
    !
    real                         :: disbk    ! local bulk dissipation
    real, dimension(mkyc,mkxc)   :: dwdx     ! x-derivative of Wigner distribution
    real, dimension(mkyc,mkxc)   :: dwdy     ! y-derivative of Wigner distribution
    real, dimension(mkyc,mkxc,5) :: W        ! the Wigner distribution in wave number space
    !
    logical                      :: stpnow   ! indicate whether program must be terminated or not
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'QCSOURCE')
    !
    if ( IGEN /= 4 ) goto 10
    !
    ! set relevant matrix elements to 0
    !
    if ( OPTG /= 5 ) then
       imatra = 0.
       imatda = 0.
    endif
    !
    ! set all dissipation coeff to 0
    !
    dissc0(1:MDC,1:MSC,1:MDISP) = 0.
    dissc1(1:MDC,1:MSC,1:MDISP) = 0.
    !
    ! set all generation coeff to 0
    !
    genc0(1:MDC,1:MSC,1:MGENR) = 0.
    genc1(1:MDC,1:MSC,1:MGENR) = 0.
    !
    ! set all redistribution coeff to 0
    !
    redc0(1:MDC,1:MSC,1:MREDS) = 0.
    redc1(1:MDC,1:MSC,1:MREDS) = 0.
    !
    ! compute surf breaking
    !
!TIMG    call SWTSTA(131)
    if ( ISURF > 0 ) then
       !
       ! calculate the quasi-homogeneous surf breaking in every sweep for the
       ! first iteration and also compute the bulk dissipation
       !
       call SSURF ( etot  , hm    , qb    , smebrk, kteta , kmespc, spcsig, ac2   ,  &
                    imatra, imatda, idcmin, idcmax, plwbrk,                          &
                    isstop, dissc0, dissc1, disbk , iter  )
       !
       ! calculate the QC surf breaking for all sweeps together from the
       ! second iteration onwards
       !
       if ( swpdir == 1 .or.                                      &
          ( swpdir == 2 .and. ix == 1 ) .or.                      &
          ( swpdir == 3 .and. iy == 1 ) .or.                      &
          ( swpdir == 4 .and. (ix == MXC .and. iy == 1 ) ) ) then
          !
          ! store bulk dissipation at current iteration
          !
          disbk1(KCGRD(1)) = disbk
          !
          ! after first iteration spatial distribution of bulk dissipation is available
          ! for computing discrete Fourier transforms
          !
          if ( iter > 1 ) then
             !
             call SWQCSURF ( memqcb, ac2, dep2, cfd, wfd, wsavd, kwave, cgo, spcdir, spcsig )
             !
          endif
          !
       endif
       !
    endif
!TIMG    call SWTSTO(131)
    !
 10 continue
    !
    ! compute quasi-coherent interactions due to medium (depth, current)
    !
!TIMG    call SWTSTA(146)
    if ( IQCM > 0 ) then
       !
       ! calculate the QC scattering for all sweeps together
       !
       if ( swpdir == 1 .or.                                      &
          ( swpdir == 2 .and. ix == 1 ) .or.                      &
          ( swpdir == 3 .and. iy == 1 ) .or.                      &
          ( swpdir == 4 .and. (ix == MXC .and. iy == 1 ) ) ) then
          !
          ! first, perform discrete Fourier transforms of the modulations ...
          !
          if ( IQCM == 1 ) call SWQCDFT ( sigft, cgft, dep2, kwave(1,1), cgo(1,1), cft, rft, sft, wft, wsave )
          if ( ICUR == 1 ) call SWQCUFT ( uxft , uyft, dep2, ux2       , uy2     , cft, rft, sft, wft, wsave )
          if ( stpnow() ) return
          !
          ! next, compute the Wigner distribution and its derivatives in geographical space ...
          !
          if ( OPTG /= 5 ) then
             ! structured mesh
             call SWQCWIG ( W, dwdx, dwdy, ac2, dep2, rdx, rdy, spcdir, spcsig )
          else
             ! unstructured mesh
             call SwanGradWig ( W, dwdx, dwdy, ac2, dep2, spcdir, spcsig )
          endif
          !
          ! ... and finally, compute the scatterer
          !
          call SWQCSCAT ( memqcm, W(1,1,1), dwdx, dwdy, sigft, cgft, uxft, uyft, dep2, kwave, cgo, spcdir, spcsig )
          !
       endif
       !
    endif
!TIMG    call SWTSTO(146)
    !
    ! get source term values for the bin that fall within a sweep and store in right hand vector
    !
    call FILQCM ( imatra, idcmin, idcmax, isstop, memqcm, memqcb, plqcs, plwbrk, redc0, dissc0 )
    !
end subroutine QCSOURCE
!
subroutine SWQCWIG ( W, dwdx, dwdy, ac2, dep2, rdx, rdy, spcdir, spcsig )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!   41.90, June 2021: New subroutine
!
!   Purpose
!
!   Computes the Wigner spectrum and its spatial derivatives on regular grid
!
!   Method
!
!   The spatial derivatives are approximated with central differences, uing 5-point stencil
!   but fall back to the first order scheme in case of dry points, boundaries, etc.
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use swcomm4
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,MCGRD), intent(in) :: ac2    ! action density at current time level
    real, dimension(MCGRD), intent(in)         :: dep2   ! water depth at current time level
    real, dimension(mkyc,mkxc), intent(out)    :: dwdx   ! x-derivative of Wigner distribution
    real, dimension(mkyc,mkxc), intent(out)    :: dwdy   ! y-derivative of Wigner distribution
    real, dimension(2), intent(in)             :: rdx    ! first component of contravariant base vector rdx(b) = a^(b)_1
    real, dimension(2), intent(in)             :: rdy    ! second component of contravariant base vector rdy(b) = a^(b)_2
    real, dimension(MDC,6), intent(in)         :: spcdir ! (*,1): spectral direction bins (radians)
                                                         ! (*,2): cosine of spectral directions
                                                         ! (*,3): sine of spectral directions
                                                         ! (*,4): cosine^2 of spectral directions
                                                         ! (*,5): cosine*sine of spectral directions
                                                         ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)           :: spcsig ! relative frequency bins
    real, dimension(mkyc,mkxc,5), intent(out)  :: W      ! the Wigner distribution in wave number space
!
!   Local variables
!
    integer       :: ic       ! loop counter over stencil
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: ik       ! index of frequency bin
    integer       :: ikx      ! loop counter over x-components of wave number space
    integer       :: iky      ! loop counter over y-components of wave number space
    integer       :: indx     ! index of grid point
    integer       :: jk       ! index of direction bin
    !
    real          :: cg       ! group velocity in wave number space
    real          :: ctb      ! contribution term
    real          :: dir      ! spectral directions in wave number space
    real          :: dp       ! local depth
    real          :: ikb      ! broken coordinate for given point in frequency-direction space
    real          :: jkb      ! broken coordinate for given point in frequency-direction space
    real          :: k        ! absolute wave number in wave number space
    real          :: sig      ! angular frequencies in wave number space
    real          :: wi1      ! first weight factor for distance in sigma-direction
    real          :: wi2      ! second weight factor for distance in sigma-direction
    real          :: wj1      ! first weight factor for distance in theta-direction
    real          :: wj2      ! second weight factor for distance in theta-direction
    !
    logical       :: outsid   ! indicates if interpolated value is/is not in grid
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCWIG')
    !
    ! first, transform action density to Wigner distribution in wave number space ...
    !
    do ic = 1, 5     ! central differences with 5-point stencil
       !
       indx = KCGRD(ic)
       !
       dp = dep2(indx)
       !
       if ( dp > DEPMIN .and. indx > 1 ) then
          !
          do ikx = 1, mkxc
             do iky = 1, mkyc
                !
                outsid = .false.
                !
                ! calculate absolute wave number, frequency, direction and group velocity in wave number space
                !
                k   = sqrt( kx(ikx)**2 + ky(iky)**2 )
                if ( k > 0. ) then
                   sig = sqrt( k * GRAV * tanh(k*dp) )
                   dir = atan2 ( ky(iky), kx(ikx) )
                   cg  = ( sig / k ) * 0.5 * ( 1. + 2.*k*dp / ( sinh(2.*k*dp) ) )
                else
                   k   = 1.
                   sig = spcsig(1)
                   dir = spcdir(1,1)
                   cg  = 0.
                   outsid = .true.
                endif
                !
                ! interpolate action density to given frequency and direction
                !
                ikb = log( sig/spcsig(1) ) / FRINTF     ! frequencies are logarithmically distributed
                !
                if ( ikb < 0. ) then
                   outsid = .true.
                else if ( ikb > real(MSC-1) ) then
                   outsid = .true.
                else if ( .not. ikb /= real(MSC-1) ) then
                   ik  = MSC - 1
                   wi2 = 1.
                else
                   ik  = int(ikb)
                   wi2 = ikb - real(ik)
                   ik  = ik + 1
                endif
                !
                jkb = ( dir - spcdir(1,1) ) / DDIR
                !
                if ( jkb < 0. ) then
                   outsid = .true.
                else if ( jkb > real(MDC-1) ) then
                   outsid = .true.
                else if ( .not. jkb /= real(MDC-1) ) then
                   jk  = MDC - 1
                   wj2 = 1.
                else
                   jk  = int(jkb)
                   wj2 = jkb - real(jk)
                   jk  = jk + 1
                endif
                !
                if ( outsid ) then
                   ctb = 0.
                else
                   wi1 = 1.- wi2
                   wj1 = 1.- wj2
                   ctb = wi1*wj1*ac2(jk,ik,indx) + wi1*wj2*ac2(jk+1,ik,indx) + wi2*wj1*ac2(jk,ik+1,indx) + wi2*wj2*ac2(jk+1,ik+1,indx)
                endif
                !
                ! multiply with the Jacobian
                !
                W(iky,ikx,ic) = ctb * cg / k
                !
             enddo
          enddo
          !
       else
          !
          W(:,:,ic) = 0.
          !
       endif
       !
    enddo
    !
    ! ... then compute the derivatives
    !
    if (PROPSL == 2 ) then
       !
       do ikx = 1, mkxc
          do iky = 1, mkyc
             !
             if ( W(iky,ikx,1) /= 0. ) then
                !
                dwdx(iky,ikx) = 0.5*rdx(1) * ( W(iky,ikx,5) - W(iky,ikx,2) ) + 0.5*rdx(2) * ( W(iky,ikx,4) - W(iky,ikx,3) )
                !
                dwdy(iky,ikx) = 0.5*rdy(1) * ( W(iky,ikx,5) - W(iky,ikx,2) ) + 0.5*rdy(2) * ( W(iky,ikx,4) - W(iky,ikx,3) )
                !
             else
                !
                dwdx(iky,ikx) = 0.
                !
                dwdy(iky,ikx) = 0.
                !
             endif
             !
          enddo
       enddo
       !
    elseif (PROPSL == 1 ) then
       !
       do ikx = 1, mkxc
          do iky = 1, mkyc
             !
             if ( W(iky,ikx,1) /= 0. ) then
                !
                dwdx(iky,ikx) = (rdx(1)+rdx(2))*W(iky,ikx,1) - rdx(1)*W(iky,ikx,2) - rdx(2)*W(iky,ikx,3)
                !
                dwdy(iky,ikx) = (rdy(1)+rdy(2))*W(iky,ikx,1) - rdy(1)*W(iky,ikx,2) - rdy(2)*W(iky,ikx,3)
                !
             else
                !
                dwdx(iky,ikx) = 0.
                !
                dwdy(iky,ikx) = 0.
                !
             endif
             !
          enddo
       enddo
       !
    endif
    !
end subroutine SWQCWIG
!
subroutine SwanGradWig ( W, dwdx, dwdy, ac2, dep2, spcdir, spcsig )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!   41.90, November 2021: New subroutine
!
!   Purpose
!
!   Computes the Wigner spectrum and its spatial derivatives on unstructured mesh
!
!   Method
!
!   The spatial derivatives are approximated by means of the Green-Gauss gradient formula
!   with the assumption of a constant gradient over the controle volume (centroid dual)
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real, dimension(nverts), intent(in)         :: dep2   ! water depth at current time level
    real, dimension(mkyc,mkxc), intent(out)     :: dwdx   ! x-derivative of Wigner distribution
    real, dimension(mkyc,mkxc), intent(out)     :: dwdy   ! y-derivative of Wigner distribution
    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
                                                          ! (*,2): cosine of spectral directions
                                                          ! (*,3): sine of spectral directions
                                                          ! (*,4): cosine^2 of spectral directions
                                                          ! (*,5): cosine*sine of spectral directions
                                                          ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
    real, dimension(mkyc,mkxc,1), intent(out)   :: W      ! the Wigner distribution in wave number space
!
!   Local variables
!
    integer                      :: icell    ! index of present cell
    integer, save                :: ient = 0 ! number of entries in this subroutine
    integer                      :: ik       ! index of frequency bin
    integer                      :: ikx      ! loop counter over x-components of wave number space
    integer                      :: iky      ! loop counter over y-components of wave number space
    integer                      :: ivert    ! index of present vertex
    integer                      :: j        ! loop counter
    integer                      :: jc       ! loop counter over cells
    integer                      :: jcell    ! index of next cell
    integer                      :: jk       ! index of direction bin
    !
    integer, dimension(3)        :: v        ! vertices in present cell
    !
    real*8                       :: area     ! twices the area of centroid dual around present vertex
    real                         :: cg       ! group velocity in wave number space
    real                         :: ctb      ! contribution term
    real                         :: dir      ! spectral directions in wave number space
    real, dimension(3)           :: dloc     ! local depth at vertices
    real                         :: ikb      ! broken coordinate for given point in frequency-direction space
    real                         :: jkb      ! broken coordinate for given point in frequency-direction space
    real                         :: k        ! absolute wave number in wave number space
    real                         :: sig      ! angular frequencies in wave number space
    real, dimension(mkyc,mkxc)   :: w0       ! Wigner spectrum in centroid of present cell
    real, dimension(mkyc,mkxc)   :: w1       ! Wigner spectrum in centroid of next cell
    real                         :: wi1      ! first weight factor for distance in sigma-direction
    real                         :: wi2      ! second weight factor for distance in sigma-direction
    real                         :: wj1      ! first weight factor for distance in theta-direction
    real                         :: wj2      ! second weight factor for distance in theta-direction
    real, dimension(mkyc,mkxc,3) :: wloc     ! local Wigner spectrum at vertices
    real*8                       :: x0       ! x-coordinate of the centroid of present cell
    real*8                       :: x1       ! x-coordinate of the centroid of next cell
    real*8                       :: y0       ! y-coordinate of the centroid of present cell
    real*8                       :: y1       ! y-coordinate of the centroid of next cell
    !
    logical                      :: outsid   ! indicates if interpolated value is/is not in grid
    logical                      :: vp       ! true if present vertex is already considered
    !
    character(80)                :: msgstr   ! string to pass message
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
    if (ltrace) call strace (ient,'SwanGradWig')
    !
    W    = 0.
    dwdx = 0.
    dwdy = 0.
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ivert = vs(1)
    !
    if ( vert(ivert)%atti(VMARKER) == 1 ) return    ! boundary vertex
    !
    area = 0.
    vp   = .false.
    !
    ! loop over cells around considered vertex
    !
    do jc = 1, vert(ivert)%noc
       !
       ! get present cell and its vertices
       !
       icell = vert(ivert)%cell(jc)%atti(CELLID)
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       dloc(1) = dep2(v(1))
       dloc(2) = dep2(v(2))
       dloc(3) = dep2(v(3))
       !
       if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
       !
       ! transform action density to Wigner distribution in wave number space for present cell
       !
       do j = 1, 3
          !
          if ( v(j) == ivert .and. vp ) then
             wloc(:,:,j) = W(:,:,1)
             cycle
          endif
          !
          do ikx = 1, mkxc
             do iky = 1, mkyc
                !
                outsid = .false.
                !
                ! calculate absolute wave number, frequency, direction and group velocity in wave number space
                !
                k   = sqrt( kx(ikx)**2 + ky(iky)**2 )
                if ( k > 0. ) then
                   sig = sqrt( k * GRAV * tanh(k*dloc(j)) )
                   dir = atan2 ( ky(iky), kx(ikx) )
                   cg  = ( sig / k ) * 0.5 * ( 1. + 2.*k*dloc(j) / ( sinh(2.*k*dloc(j)) ) )
                else
                   k   = 1.
                   sig = spcsig(1)
                   dir = spcdir(1,1)
                   cg  = 0.
                   outsid = .true.
                endif
                !
                ! interpolate action density to given frequency and direction
                !
                ikb = log( sig/spcsig(1) ) / FRINTF     ! frequencies are logarithmically distributed
                !
                if ( ikb < 0. ) then
                   outsid = .true.
                else if ( ikb > real(MSC-1) ) then
                   outsid = .true.
                else if ( .not. ikb /= real(MSC-1) ) then
                   ik  = MSC - 1
                   wi2 = 1.
                else
                   ik  = int(ikb)
                   wi2 = ikb - real(ik)
                   ik  = ik + 1
                endif
                !
                jkb = ( dir - spcdir(1,1) ) / DDIR
                !
                if ( jkb < 0. ) then
                   outsid = .true.
                else if ( jkb > real(MDC-1) ) then
                   outsid = .true.
                else if ( .not. jkb /= real(MDC-1) ) then
                   jk  = MDC - 1
                   wj2 = 1.
                else
                   jk  = int(jkb)
                   wj2 = jkb - real(jk)
                   jk  = jk + 1
                endif
                !
                if ( outsid ) then
                   ctb = 0.
                else
                   wi1 = 1.- wi2
                   wj1 = 1.- wj2
                   ctb = wi1*wj1*ac2(jk,ik,v(j)) + wi1*wj2*ac2(jk+1,ik,v(j)) + wi2*wj1*ac2(jk,ik+1,v(j)) + wi2*wj2*ac2(jk+1,ik+1,v(j))
                endif
                !
                ! multiply with the Jacobian
                !
                wloc(iky,ikx,j) = ctb * cg / k
                !
             enddo
          enddo
          !
          ! store Wigner spectrum in present vertex
          !
          if ( v(j) == ivert .and. .not.vp ) then
             W(:,:,1) = wloc(:,:,j)
             vp       = .true.
          endif
          !
       enddo
       !
       ! determine centroid of present cell
       !
       x0 = cell(icell)%attr(CELLCX)
       y0 = cell(icell)%attr(CELLCY)
       !
       ! determine Wigner spectrum in centroid in present cell
       !
       w0(:,:) = ( wloc(:,:,1) + wloc(:,:,2) + wloc(:,:,3) )/ 3.
       !
       ! get next cell in counterclockwise direction
       !
       jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
       !
       v(1) = cell(jcell)%atti(CELLV1)
       v(2) = cell(jcell)%atti(CELLV2)
       v(3) = cell(jcell)%atti(CELLV3)
       !
       dloc(1) = dep2(v(1))
       dloc(2) = dep2(v(2))
       dloc(3) = dep2(v(3))
       !
       if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
       !
       ! transform action density to Wigner distribution in wave number space for next cell
       !
       do j = 1, 3
          !
          if ( v(j) == ivert .and. vp ) then
             wloc(:,:,j) = W(:,:,1)
             cycle
          endif
          !
          do ikx = 1, mkxc
             do iky = 1, mkyc
                !
                outsid = .false.
                !
                ! calculate absolute wave number, frequency, direction and group velocity in wave number space
                !
                k   = sqrt( kx(ikx)**2 + ky(iky)**2 )
                if ( k > 0. ) then
                   sig = sqrt( k * GRAV * tanh(k*dloc(j)) )
                   dir = atan2 ( ky(iky), kx(ikx) )
                   cg  = ( sig / k ) * 0.5 * ( 1. + 2.*k*dloc(j) / ( sinh(2.*k*dloc(j)) ) )
                else
                   k   = 1.
                   sig = spcsig(1)
                   dir = spcdir(1,1)
                   cg  = 0.
                   outsid = .true.
                endif
                !
                ! interpolate action density to given frequency and direction
                !
                ikb = log( sig/spcsig(1) ) / FRINTF     ! frequencies are logarithmically distributed
                !
                if ( ikb < 0. ) then
                   outsid = .true.
                else if ( ikb > real(MSC-1) ) then
                   outsid = .true.
                else if ( .not. ikb /= real(MSC-1) ) then
                   ik  = MSC - 1
                   wi2 = 1.
                else
                   ik  = int(ikb)
                   wi2 = ikb - real(ik)
                   ik  = ik + 1
                endif
                !
                jkb = ( dir - spcdir(1,1) ) / DDIR
                !
                if ( jkb < 0. ) then
                   outsid = .true.
                else if ( jkb > real(MDC-1) ) then
                   outsid = .true.
                else if ( .not. jkb /= real(MDC-1) ) then
                   jk  = MDC - 1
                   wj2 = 1.
                else
                   jk  = int(jkb)
                   wj2 = jkb - real(jk)
                   jk  = jk + 1
                endif
                !
                if ( outsid ) then
                   ctb = 0.
                else
                   wi1 = 1.- wi2
                   wj1 = 1.- wj2
                   ctb = wi1*wj1*ac2(jk,ik,v(j)) + wi1*wj2*ac2(jk+1,ik,v(j)) + wi2*wj1*ac2(jk,ik+1,v(j)) + wi2*wj2*ac2(jk+1,ik+1,v(j))
                endif
                !
                ! multiply with the Jacobian
                !
                wloc(iky,ikx,j) = ctb * cg / k
                !
             enddo
          enddo
          !
       enddo
       !
       ! determine centroid of next cell
       !
       x1 = cell(jcell)%attr(CELLCX)
       y1 = cell(jcell)%attr(CELLCY)
       !
       ! determine Wigner spectrum in centroid of next cell
       !
       w1(:,:) = ( wloc(:,:,1) + wloc(:,:,2) + wloc(:,:,3) )/ 3.
       !
       ! compute contribution to area of centroid dual
       !
       area = area + x0*y1 - x1*y0
       !
       ! compute contribution to x-gradient of Wigner spectrum
       !
       dwdx(:,:) = dwdx(:,:) + ( w0(:,:) + w1(:,:) ) * real( y1 - y0 )
       !
       ! compute contribution to y-gradient of Wigner spectrum
       !
       dwdy(:,:) = dwdy(:,:) + ( w0(:,:) + w1(:,:) ) * real( x1 - x0 )
       !
    enddo
    !
    ! if area is non-positive, give error and go to next vertex
    !
    if ( .not. area > 0. ) then
       write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
       call msgerr( 2, trim(msgstr) )
       return
    endif
    !
    dwdx(:,:) =  dwdx(:,:)/real(area)
    dwdy(:,:) = -dwdy(:,:)/real(area)
    !
    do ikx = 1, mkxc
       do iky = 1, mkyc
          if ( .not. W(iky,ikx,1) /= 0. ) then
             dwdx(iky,ikx) = 0.
             dwdy(iky,ikx) = 0.
          endif
       enddo
    enddo
    !
    return
    !
 10 W    = 0.
    dwdx = 0.
    dwdy = 0.
    !
end subroutine SwanGradWig
!
subroutine SWQCSCAT ( memqcm, W, dwdx, dwdy, sigft, cgft, uxft, uyft, dep2, kwave, cgo, spcdir, spcsig )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!   41.90, June 2021: New subroutine
!
!   Purpose
!
!   Computation of the source term due to QC scattering for all sweeps together
!
!   Method
!
!   See e.g. the JPO paper of Smit et al (2015), Eq. (15) and the JFM paper of Akrish et al (2020), Eq. (2.21)
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
!
    implicit none
!
!   Argument variables
!
    real   , dimension(ncoz,ncoz,MSC), intent(in)  :: cgft   ! Fourier-transformed modulation of group velocity
    real   , dimension(MSC,MICMAX)   , intent(in)  :: cgo    ! group velocity
    real   , dimension(MCGRD)        , intent(in)  :: dep2   ! water depth at current time level
    real   , dimension(mkyc,mkxc)    , intent(in)  :: dwdx   ! x-derivative of Wigner distribution
    real   , dimension(mkyc,mkxc)    , intent(in)  :: dwdy   ! y-derivative of Wigner distribution
    real   , dimension(MSC,MICMAX)   , intent(in)  :: kwave  ! wave number
    real   , dimension(MDC,MSC,MCGRD), intent(out) :: memqcm ! auxiliary array to store results of QC scattering in full spectral space
    real   , dimension(ncoz,ncoz,MSC), intent(in)  :: sigft  ! Fourier-transformed modulation of intrinsic frequency
    real   , dimension(MDC,6)        , intent(in)  :: spcdir ! (*,1): spectral direction bins (radians)
                                                             ! (*,2): cosine of spectral directions
                                                             ! (*,3): sine of spectral directions
                                                             ! (*,4): cosine^2 of spectral directions
                                                             ! (*,5): cosine*sine of spectral directions
                                                             ! (*,6): sine^2 of spectral directions
    real   , dimension(MSC)          , intent(in) :: spcsig  ! relative frequency bins
    complex, dimension(ncoz,ncoz)    , intent(in) :: uxft    ! u-component of Fourier-transformed modulation of ambient current
    complex, dimension(ncoz,ncoz)    , intent(in) :: uyft    ! v-component of Fourier-transformed modulation of ambient current
    real   , dimension(mkyc,mkxc)    , intent(in) :: W       ! the Wigner distribution in wave number space
!
!   Local variables
!
    integer                    :: id       ! loop counter over direction bins
    integer, save              :: ient = 0 ! number of entries in this subroutine
    integer                    :: ik       ! loop/offset index of kx
    integer                    :: ikm      ! offset index for interpolation
    integer                    :: ikx      ! loop counter over x-components of wave number space
    integer                    :: iky      ! loop counter over y-components of wave number space
    integer                    :: is       ! loop counter over frequency bins
    integer                    :: ish      ! ish+1 corresponds to the centre of coherent scattering zone
    integer                    :: jk       ! loop/offset index of ky
    !
    real                       :: cax      ! wave transport velocity in x-direction
    real                       :: cay      ! wave transport velocity in y-direction
    real                       :: ctb      ! contribution term
    real                       :: cg       ! local value of Fourier-transformed modulation of group velocity
    real                       :: dkx      ! resolution of wave number space in x-direction
    real                       :: dky      ! resolution of wave number space in y-direction
    real                       :: dp       ! local depth
    real                       :: ikb      ! broken coordinate for given point in wave number space
    real                       :: jkb      ! broken coordinate for given point in wave number space
    real                       :: k        ! local value of absolute wave number
    real                       :: kxl      ! x-component of local wave number
    real                       :: kyl      ! y-component of local wave number
    real                       :: sig      ! relative frequency / local value of Fourier-transformed modulation of intrinsic frequency
    real                       :: wf1      ! first interpolation weight factor
    real                       :: wf2      ! second interpolation weight factor
    real                       :: wi1      ! first weight factor for distance in kx-direction
    real                       :: wi2      ! second weight factor for distance in kx-direction
    real                       :: wj1      ! first weight factor for distance in ky-direction
    real                       :: wj2      ! second weight factor for distance in ky-direction
    !
    complex                    :: u        ! local u-component of Fourier-transformed modulation
    complex                    :: v        ! local v-component of Fourier-transformed modulation
    !
    real, dimension(mkyc,mkxc) :: arr      ! array to store data temporarily
    !
    logical                    :: outsid   ! indicates if interpolated value is/is not in grid
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCSCAT')
    !
    memqcm(:,:,KCGRD(1)) = 0.
    !
    dp = dep2(KCGRD(1))
    if ( .not. dp > DEPMIN ) return
    !
    ish = ncoz/2
    !
    do ikx = 1, mkxc
       !
       kxl = kx(ikx)
       !
       do iky = 1, mkyc
          !
          kyl = ky(iky)
          !
          k = sqrt( kxl*kxl + kyl*kyl )
          !
          if ( k > 0. ) then
             sig = sqrt( k * GRAV * tanh(k*dp) )
          else
             k   = 1.
             sig = spcsig(1)
          endif
          !
          ikb = log( sig/spcsig(1) ) / FRINTF     ! frequencies are logarithmically distributed
          !
          if ( ikb < 0. ) then
             ikm = 1
             wf2 = 0.
          else if ( .not. ikb < real(MSC-1) ) then
             ikm = MSC - 1
             wf2 = 1.
          else
             ikm = int(ikb)
             wf2 = ikb - real(ikm)
             ikm = ikm + 1
          endif
          !
          wf1 = 1. - wf2
          !
          cg = wf1 * cgft(ish+1,ish+1,ikm) + wf2 * cgft(ish+1,ish+1,ikm+1)
          !
          u  = uxft(ish+1,ish+1)
          v  = uyft(ish+1,ish+1)
          !
          cax = cg * kxl / k + real(u)
          cay = cg * kyl / k + real(v)
          !
          ctb = cax * dwdx(iky,ikx) + cay * dwdy(iky,ikx)
          !
          ! sum in the negative direction of kx and ky
          !
          do ik = 0, min(ikx-1,ish-1)
             do jk = 0, min(iky-1,ish-1)
                !
                sig = wf1 * sigft(ish+1+ik,ish+1+jk,ikm) + wf2 * sigft(ish+1+ik,ish+1+jk,ikm+1)
                cg  = wf1 * cgft (ish+1+ik,ish+1+jk,ikm) + wf2 * cgft (ish+1+ik,ish+1+jk,ikm+1)
                !
                u   = uxft(ish+1+ik,ish+1+jk)
                v   = uyft(ish+1+ik,ish+1+jk)
                !
                cax = cg * kxl / k + real(u)
                cay = cg * kyl / k + real(v)
                !
                ctb = ctb + 2. * ( sig + kxl*imag(u) + kyl*imag(v) ) * W(iky-jk,ikx-ik) - cax * dwdx(iky-jk,ikx-ik) - cay * dwdy(iky-jk,ikx-ik)
                !
             enddo
          enddo
          !
          ! sum in the positive direction of kx and ky
          !
          do ik = 0, min(mkxc-ikx,ish-1)
             do jk = 0, min(mkyc-iky,ish-1)
                !
                sig = wf1 * sigft(ish+1+ik,ish+1+jk,ikm) + wf2 * sigft(ish+1+ik,ish+1+jk,ikm+1)
                cg  = wf1 * cgft (ish+1+ik,ish+1+jk,ikm) + wf2 * cgft (ish+1+ik,ish+1+jk,ikm+1)
                !
                u   = uxft(ish+1+ik,ish+1+jk)
                v   = uyft(ish+1+ik,ish+1+jk)
                !
                cax = cg * kxl / k + real(u)
                cay = cg * kyl / k + real(v)
                !
                ctb = ctb - 2. * ( sig + kxl*imag(u) + kyl*imag(v) ) * W(iky+jk,ikx+ik) - cax * dwdx(iky+jk,ikx+ik) - cay * dwdy(iky+jk,ikx+ik)
                !
             enddo
          enddo
          !
          ! sum in the negative direction of kx and in the positive direction of ky
          !
          do ik = 1, min(ikx-1,ish-1)
             do jk = 1, min(mkyc-iky,ish-1)
                !
                sig = wf1 * sigft(ish+1+ik,ish+1-jk,ikm) + wf2 * sigft(ish+1+ik,ish+1-jk,ikm+1)
                cg  = wf1 * cgft (ish+1+ik,ish+1-jk,ikm) + wf2 * cgft (ish+1+ik,ish+1-jk,ikm+1)
                !
                u   = uxft(ish+1+ik,ish+1-jk)
                v   = uyft(ish+1+ik,ish+1-jk)
                !
                cax = cg * kxl / k + real(u)
                cay = cg * kyl / k + real(v)
                !
                ctb = ctb + 2. * ( sig + kxl*imag(u) + kyl*imag(v) ) * W(iky+jk,ikx-ik) - cax * dwdx(iky+jk,ikx-ik) - cay * dwdy(iky+jk,ikx-ik)
                !
             enddo
          enddo
          !
          ! sum in the positive direction of kx and in the negative direction of ky
          !
          do ik = 1, min(mkxc-ikx,ish-1)
             do jk = 1, min(iky-1,ish-1)
                !
                sig = wf1 * sigft(ish+1+ik,ish+1-jk,ikm) + wf2 * sigft(ish+1+ik,ish+1-jk,ikm+1)
                cg  = wf1 * cgft (ish+1+ik,ish+1-jk,ikm) + wf2 * cgft (ish+1+ik,ish+1-jk,ikm+1)
                !
                u   = uxft(ish+1+ik,ish+1-jk)
                v   = uyft(ish+1+ik,ish+1-jk)
                !
                cax = cg * kxl / k + real(u)
                cay = cg * kyl / k + real(v)
                !
                ctb = ctb - 2. * ( sig + kxl*imag(u) + kyl*imag(v) ) * W(iky-jk,ikx+ik) - cax * dwdx(iky-jk,ikx+ik) - cay * dwdy(iky-jk,ikx+ik)
                !
             enddo
          enddo
          !
          ! store the total source term
          !
          arr(iky,ikx) = ctb
          !
       enddo
       !
    enddo
    !
    ! transform the scatterer back to frequency-direction space
    !
    dkx = ( kx(mkxc) - kx(1) ) / float(mkxc-1)
    dky = ( ky(mkyc) - ky(1) ) / float(mkyc-1)
    !
    do is = 1, MSC
       do id = 1, MDC
          !
          ! for each frequency and direction compute the corresponding wave number vector ...
          !
          kxl = kwave(is,1) * spcdir(id,2)
          kyl = kwave(is,1) * spcdir(id,3)
          !
          ! ... and subsequently the corresponding interpolation factors
          !
          ikb = ( kxl - kx(1) ) / dkx
          jkb = ( kyl - ky(1) ) / dky
          !
          outsid = .false.
          !
          if ( ikb < 0. ) then
             outsid = .true.
          else if ( ikb > real(mkxc-1) ) then
             outsid = .true.
          else if ( .not. ikb /= real(mkxc-1) ) then
             ik  = mkxc - 1
             wi2 = 1.
          else
             ik  = int(ikb)
             wi2 = ikb - real(ik)
             ik  = ik + 1
          endif
          !
          if ( jkb < 0. ) then
             outsid = .true.
          else if ( jkb > real(mkyc-1) ) then
             outsid = .true.
          else if ( .not. jkb /= real(mkyc-1) ) then
             jk  = mkyc - 1
             wj2 = 1.
          else
             jk  = int(jkb)
             wj2 = jkb - real(jk)
             jk  = jk + 1
          endif
          !
          if ( outsid ) then
             ctb = 0.
          else
             wi1 = 1.- wi2
             wj1 = 1.- wj2
             ctb = wi1*wj1*arr(jk,ik) + wi1*wj2*arr(jk+1,ik) + wi2*wj1*arr(jk,ik+1) + wi2*wj2*arr(jk+1,ik+1)
          endif
          !
          ! multiply with the inverse of Jacobian
          !
          memqcm(id,is,KCGRD(1)) = ctb * kwave(is,1) / cgo(is,1)
          !
       enddo
    enddo
    !
end subroutine SWQCSCAT
!
subroutine SWQCSURF ( memqcb, ac2, dep2, cfd, wfd, wsavd, kwave, cgo, spcdir, spcsig )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.91: Pieter Smit and Marcel Zijlema
!
!   Updates
!
!   41.91, February 2022: New subroutine
!
!   Purpose
!
!   Computes the QC source term due to surf breaking
!
!   Method
!
!   See the OM paper of Smit et al (2015), Section 2.2
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,MCGRD), intent(in)         :: ac2    ! action density at current time level
    complex(kind=8), dimension(myd,mxd), intent(inout) :: cfd    ! Fourier coefficients (FFT)
    real, dimension(MSC,MICMAX), intent(in)            :: cgo    ! group velocity
    real, dimension(MCGRD), intent(in)                 :: dep2   ! water depth at current time level
    real, dimension(MSC,MICMAX), intent(in)            :: kwave  ! wave number
    real, dimension(MDC,MSC,MCGRD), intent(out)        :: memqcb ! auxiliary array to store results of QC surf breaking in full spectral space
    real, dimension(MDC,6), intent(in)                 :: spcdir ! (*,1): spectral direction bins (radians)
                                                                 ! (*,2): cosine of spectral directions
                                                                 ! (*,3): sine of spectral directions
                                                                 ! (*,4): cosine^2 of spectral directions
                                                                 ! (*,5): cosine*sine of spectral directions
                                                                 ! (*,6): sine^2 of spectral directions
    real, dimension(MSC), intent(in)                   :: spcsig ! relative frequency bins
    real(kind=8), dimension(lenwfd), intent(inout)     :: wfd    ! work array (FFT)
    real(kind=8), dimension(lensvd), intent(inout)     :: wsavd  ! work array (FFT)
!
!   Local variables
!
    integer                  :: id       ! loop counter over direction bins
    integer, save            :: ient = 0 ! number of entries in this subroutine
    integer                  :: ierr     ! error indicator
    integer                  :: ik       ! index of frequency bin
    integer                  :: is       ! loop counter over frequency bins
    integer                  :: ixd      ! loop counter over x-components of wave number/spatial lag grid
    integer                  :: ixnyq    ! point index in x-direction corresponding to Nyquist frequency
    integer                  :: ixs      ! shifted point index in x-direction
    integer                  :: iyd      ! loop counter over y-components of wave number/spatial lag grid
    integer                  :: iynyq    ! point index in y-direction corresponding to Nyquist frequency
    integer                  :: iys      ! shifted point index in y-direction
    integer                  :: jk       ! index of direction bin
    !
    real                     :: cg       ! group velocity in wave number grid
    real                     :: ctb      ! contribution term
    real                     :: dir      ! spectral directions in wave number grid
    real                     :: disbk    ! local bulk dissipation at current grid point (from previous iteration)
    real                     :: dkx      ! resolution of wave number grid in x-direction
    real                     :: dky      ! resolution of wave number grid in y-direction
    real                     :: dp       ! local depth
    real                     :: ikb      ! broken coordinate for given point in frequency-direction/wave number grid
    real                     :: jkb      ! broken coordinate for given point in frequency-direction/wave number grid
    real                     :: k        ! absolute wave number in wave number grid
    real                     :: kxl      ! x-component of local wave number
    real                     :: kyl      ! y-component of local wave number
    real                     :: sig      ! angular frequencies in wave number grid
    real                     :: varD     ! variance of QC surf breaking in wave number space
    real                     :: varW     ! variance of Wigner distribution in wave number space
    real                     :: wi1      ! first weight factor for distance in sigma/kxd-direction
    real                     :: wi2      ! second weight factor for distance in sigma/kxd-direction
    real                     :: wj1      ! first weight factor for distance in theta/kyd-direction
    real                     :: wj2      ! second weight factor for distance in theta/kyd-direction
    !
    real, dimension(myd,mxd) :: disqc    ! bulk dissipation in QC framework
    real, dimension(myd,mxd) :: W        ! the Wigner distribution
    !
    logical                  :: outsid   ! indicates if interpolated value is/is not in grid
    !
    character(80)            :: msgstr   ! string to pass message
!
!   Structure
!
!   interpolate Wigner distribution to the wave number grid associated with the spatial lag grid
!   compute correlation function by means of the inverse Fourier transform of Wigner distribution
!   compute bulk dissipation as function of spatial lag
!   compute the product of bulk dissipation and correlation function
!   compute the dissipation term owing to the Fourier transform of the product to the wave number grid
!   transform the dissipation term back to frequency-direction grid
!
!   Source text
!
    if (ltrace) call strace (ient,'SWQCSURF')
    !
    memqcb(:,:,KCGRD(1)) = 0.
    !
    dp = dep2(KCGRD(1))
    if ( .not. dp > DEPMIN ) return
    !
    ierr = 0
    !
    ! interpolate Wigner distribution to the wave number grid associated with the spatial lag grid
    !
    do ixd = 1, mxd
       do iyd = 1, myd
          !
          outsid = .false.
          !
          ! calculate absolute wave number, frequency, direction and group velocity
          !
          k = sqrt( kxd(ixd)**2 + kyd(iyd)**2 )
          if ( k > 0. ) then
             sig = sqrt( k * GRAV * tanh(k*dp) )
             dir = atan2 ( kyd(iyd), kxd(ixd) )
             cg  = ( sig / k ) * 0.5 * ( 1. + 2.*k*dp / ( sinh(2.*k*dp) ) )
          else
             k   = 1.
             sig = spcsig(1)
             dir = spcdir(1,1)
             cg  = 0.
             outsid = .true.
          endif
          !
          ! interpolate action density to given frequency and direction
          !
          ikb = log( sig/spcsig(1) ) / FRINTF     ! frequencies are logarithmically distributed
          !
          if ( ikb < 0. ) then
             outsid = .true.
          else if ( ikb > real(MSC-1) ) then
             outsid = .true.
          else if ( .not. ikb /= real(MSC-1) ) then
             ik  = MSC - 1
             wi2 = 1.
          else
             ik  = int(ikb)
             wi2 = ikb - real(ik)
             ik  = ik + 1
          endif
          !
          jkb = ( dir - spcdir(1,1) ) / DDIR
          !
          if ( jkb < 0. ) then
             outsid = .true.
          else if ( jkb > real(MDC-1) ) then
             outsid = .true.
          else if ( .not. jkb /= real(MDC-1) ) then
             jk  = MDC - 1
             wj2 = 1.
          else
             jk  = int(jkb)
             wj2 = jkb - real(jk)
             jk  = jk + 1
          endif
          !
          if ( outsid ) then
             ctb = 0.
          else
             wi1 = 1.- wi2
             wj1 = 1.- wj2
             ctb = wi1*wj1*ac2(jk,ik,KCGRD(1)) + wi1*wj2*ac2(jk+1,ik,KCGRD(1)) + wi2*wj1*ac2(jk,ik+1,KCGRD(1)) + wi2*wj2*ac2(jk+1,ik+1,KCGRD(1))
          endif
          !
          ! multiply with the Jacobian
          !
          W(iyd,ixd) = ctb * cg / k
          !
       enddo
    enddo
    !
    ! compute correlation function by means of the inverse Fourier transform of Wigner distribution
    !
    cfd = cmplx(W)
    !
    call cfft2b ( myd, myd, mxd, cfd, wsavd, lensvd, wfd, lenwfd, ierr )
    if ( ierr /= 0 ) goto 10
    !
    ! compute bulk dissipation as function of spatial lag (see Eq. 16)
    !
    call bdiss
    !
    ixnyq = mxd/2 + 1
    iynyq = myd/2 + 1
    !
    disbk = disqc(iynyq,ixnyq)  ! = disbk0(KCGRD(1))
    !
    ! compute the product of bulk dissipation and correlation function (see Eq. 23)
    !
    do ixd = 1, mxd
       !
       if ( ixd < ixnyq ) then
          !
          ixs = ixd + ixnyq - 1 + mod(mxd,2)
          !
       else
          !
          ixs = ixd - ixnyq + 1
          !
       endif
       !
       do iyd = 1, myd
           !
           if ( iyd < iynyq ) then
               !
               iys = iyd + iynyq - 1 + mod(myd,2)
               !
           else
               !
               iys = iyd - iynyq + 1
               !
           endif
           !
           cfd(iys,ixs) = disqc(iyd,ixd) * cfd(iys,ixs)
           !
       enddo
       !
    enddo
    !
    ! compute the dissipation term owing to the Fourier transform of the product to the wave number grid (see Eq. 23)
    !
    call cfft2f ( myd, myd, mxd, cfd, wsavd, lensvd, wfd, lenwfd, ierr )
    if ( ierr /= 0 ) goto 10
    !
    disqc = real(cfd)     ! note: the imaginary part of cfd must be zero
    !
    ! scale dissipation to get the correct variance, which must be negative
    !
    dkx = ( kxd(mxd) - kxd(1) ) / float(mxd-1)
    dky = ( kyd(myd) - kyd(1) ) / float(myd-1)
    !
    varD = sum(disqc) * dkx * dky
    varW = sum(W)     * dkx * dky
    !
    call varchk
    !
    ! transform the dissipation term back to frequency-direction grid
    !
    do is = 1, MSC
       do id = 1, MDC
          !
          ! for each frequency and direction compute the corresponding wave number vector ...
          !
          kxl = kwave(is,1) * spcdir(id,2)
          kyl = kwave(is,1) * spcdir(id,3)
          !
          ! ... and subsequently the corresponding interpolation factors
          !
          ikb = ( kxl - kxd(1) ) / dkx
          jkb = ( kyl - kyd(1) ) / dky
          !
          outsid = .false.
          !
          if ( ikb < 0. ) then
             outsid = .true.
          else if ( ikb > real(mxd-1) ) then
             outsid = .true.
          else if ( .not. ikb /= real(mxd-1) ) then
             ik  = mxd - 1
             wi2 = 1.
          else
             ik  = int(ikb)
             wi2 = ikb - real(ik)
             ik  = ik + 1
          endif
          !
          if ( jkb < 0. ) then
             outsid = .true.
          else if ( jkb > real(myd-1) ) then
             outsid = .true.
          else if ( .not. jkb /= real(myd-1) ) then
             jk  = myd - 1
             wj2 = 1.
          else
             jk  = int(jkb)
             wj2 = jkb - real(jk)
             jk  = jk + 1
          endif
          !
          if ( outsid ) then
             ctb = 0.
          else
             wi1 = 1.- wi2
             wj1 = 1.- wj2
             ctb = wi1*wj1*disqc(jk,ik) + wi1*wj2*disqc(jk+1,ik) + wi2*wj1*disqc(jk,ik+1) + wi2*wj2*disqc(jk+1,ik+1)
          endif
          !
          ! multiply with the inverse of Jacobian
          !
          memqcb(id,is,KCGRD(1)) = ctb * kwave(is,1) / cgo(is,1)
          !
       enddo
    enddo
    !
 10 if ( ierr /= 0 ) then
       write (msgstr, '(a,i6)') 'something went wrong with the Fourier transform - return code is ',ierr
       call msgerr ( 4, trim(msgstr) )
    endif
    !
    contains
    !
    subroutine bdiss
    !
    ! note: this subroutine makes use of geographical grid (regular, curvilinear or unstructured) and associated data
    !
    use swcomm2
    use m_genarr
    use m_parall
    use SwanGriddata
    use SwanCompdata
    !
    integer, dimension(6)              :: itmp     ! temporary stored integers for swapping
    integer                            :: ixg      ! index of computational grid point in x-direction
    integer                            :: iyg      ! index of computational grid point in y-direction
    integer                            :: j        ! counter
    !
    real, dimension(4)                 :: rtmp     ! temporary stored reals for swapping
    real                               :: x        ! x-coordinate of computational point offset by spatial lag
    real                               :: xcl      ! local broken x-coordinate
    real                               :: xp       ! x-coordinate of computational grid point
    real                               :: y        ! y-coordinate of computational point offset by spatial lag
    real                               :: ycl      ! local broken y-coordinate
    real                               :: yp       ! y-coordinate of computational grid point
    !
    real, dimension(:,:), allocatable  :: dneg     ! negative contribution to bulk dissipation
    real, dimension(:,:), allocatable  :: dpos     ! positive contribution to bulk dissipation
    real, dimension(:)  , allocatable  :: xc       ! broken x-coordinate
    real, dimension(:)  , allocatable  :: yc       ! broken y-coordinate
    !
    ! allocation of temporary arrays
    !
    allocate ( xc( mxd*myd ) )
    allocate ( yc( mxd*myd ) )
    !
    allocate ( dneg( myd,mxd ) )
    allocate ( dpos( myd,mxd ) )
    !
    ! first, compute broken coordinates for the offset points to be interpolated ...
    !
    if ( optg == 1 ) then
       !
       ! regular grid
       !
       xp = xcgrid(IXCGRD(1),IYCGRD(1))
       yp = ycgrid(IXCGRD(1),IYCGRD(1))
       !
       j = 0
       !
       do ixd = 1, mxd
          !
          x = xp + 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp + 0.5 * ypd(iyd)
             !
             j = j + 1
             !
             xc(j) = (  ( x - XPC ) * COSPC + ( y - YPC ) * SINPC ) / DX
             yc(j) = ( -( x - XPC ) * SINPC + ( y - YPC ) * COSPC ) / DY
             !
          enddo
          !
       enddo
       !
    else if ( optg == 3 ) then
       !
       ! curvilinear grid
       !
       xp = xcgrid(IXCGRD(1),IYCGRD(1))
       yp = ycgrid(IXCGRD(1),IYCGRD(1))
       !
       itmp(1) = MXC
       itmp(2) = MYC
       itmp(3) = MCGRD
       itmp(4) = NGRBND
       itmp(5) = MXF
       itmp(6) = MYF
       !
       rtmp(1) = XCLMIN
       rtmp(2) = XCLMAX
       rtmp(3) = YCLMIN
       rtmp(4) = YCLMAX
       !
       MXC    = MXCGL
       MYC    = MYCGL
       MCGRD  = MCGRDGL
       NGRBND = NGRBGL
       MXF    = 1
       MYF    = 1
       !
       XCLMIN = XCGMIN
       XCLMAX = XCGMAX
       YCLMIN = YCGMIN
       YCLMAX = YCGMAX
       !
       j = 0
       !
       do ixd = 1, mxd
          !
          x = xp + 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp + 0.5 * ypd(iyd)
             !
             call CVMESH ( x, y, xcl, ycl, KGRPGL, XGRDGL ,YGRDGL, KGRBGL )
             !
             j = j + 1
             !
             xc(j) = xcl
             yc(j) = ycl
             !
          enddo
          !
       enddo
       !
       MXC    = itmp(1)
       MYC    = itmp(2)
       MCGRD  = itmp(3)
       NGRBND = itmp(4)
       MXF    = itmp(5)
       MYF    = itmp(6)
       !
       XCLMIN = rtmp(1)
       XCLMAX = rtmp(2)
       YCLMIN = rtmp(3)
       YCLMAX = rtmp(4)
       !
    endif
    !
    ! ... then carry out the interpolation for positive contribution to bulk dissipation
    !
    if ( optg /= 5 ) then
       !
       j = 0
       !
       do ixd = 1, mxd
          do iyd = 1, myd
             !
             j = j + 1
             !
             outsid = .false.
             !
             if ( xc(j) <= 0. ) then
                ixg = 1
                wi2 = 0.
                if ( xc(j) < -0.1 ) outsid = .true.
             else if ( xc(j) >= float(MXCGL-1) ) then
                ixg = MXCGL - 1
                wi2 = 1.
                if ( xc(j) > float(MXCGL)-0.9 ) outsid = .true.
             else
                ixg = int(xc(j))
                wi2 = xc(j) - real(ixg)
                ixg = ixg + 1
             endif
             !
             if ( yc(j) <= 0. ) then
                iyg = 1
                wj2 = 0.
                if ( yc(j) < -0.1 ) outsid = .true.
             else if ( yc(j) >= float(MYCGL-1) ) then
                iyg = MYCGL - 1
                wj2 = 1.
                if ( yc(j) > float(MYCGL)-0.9 ) outsid = .true.
             else
                iyg = int(yc(j))
                wj2 = yc(j) - real(iyg)
                iyg = iyg + 1
             endif
             !
             wi1 = 1.- wi2
             wj1 = 1.- wj2
             !
             if ( outsid ) then
                !
                dpos(iyd,ixd) = 0.
                !
             else
                !
                dpos(iyd,ixd) = wi1*wj1*disbk0(KGRPGL(ixg,iyg)) + wi1*wj2*disbk0(KGRPGL(ixg,iyg+1)) + wi2*wj1*disbk0(KGRPGL(ixg+1,iyg)) + wi2*wj2*disbk0(KGRPGL(ixg+1,iyg+1))
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       ! unstructured mesh
       !
       xp = xcugrd(vs(1))
       yp = ycugrd(vs(1))
       !
       do ixd = 1, mxd
          !
          x = xp + 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp + 0.5 * ypd(iyd)
             !
             if ( x < XCGMIN .or. x > XCGMAX .or. y < YCGMIN .or. y > YCGMAX ) then
                !
                dpos(iyd,ixd) = 0.
                !
             else
                !
                call SwanInterpolatePoint( dpos(iyd,ixd), x, y, disbk0, 0. )
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! repeat for the negative part of bulk dissipation ...
    !
    if ( optg == 1 ) then
       !
       ! regular grid
       !
       xp = xcgrid(IXCGRD(1),IYCGRD(1))
       yp = ycgrid(IXCGRD(1),IYCGRD(1))
       !
       j = 0
       !
       do ixd = 1, mxd
          !
          x = xp - 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp - 0.5 * ypd(iyd)
             !
             j = j + 1
             !
             xc(j) = (  ( x - XPC ) * COSPC + ( y - YPC ) * SINPC ) / DX
             yc(j) = ( -( x - XPC ) * SINPC + ( y - YPC ) * COSPC ) / DY
             !
          enddo
          !
       enddo
       !
    else if ( optg == 3 ) then
       !
       ! curvilinear grid
       !
       xp = xcgrid(IXCGRD(1),IYCGRD(1))
       yp = ycgrid(IXCGRD(1),IYCGRD(1))
       !
       itmp(1) = MXC
       itmp(2) = MYC
       itmp(3) = MCGRD
       itmp(4) = NGRBND
       itmp(5) = MXF
       itmp(6) = MYF
       !
       rtmp(1) = XCLMIN
       rtmp(2) = XCLMAX
       rtmp(3) = YCLMIN
       rtmp(4) = YCLMAX
       !
       MXC    = MXCGL
       MYC    = MYCGL
       MCGRD  = MCGRDGL
       NGRBND = NGRBGL
       MXF    = 1
       MYF    = 1
       !
       XCLMIN = XCGMIN
       XCLMAX = XCGMAX
       YCLMIN = YCGMIN
       YCLMAX = YCGMAX
       !
       j = 0
       !
       do ixd = 1, mxd
          !
          x = xp - 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp - 0.5 * ypd(iyd)
             !
             call CVMESH ( x, y, xcl, ycl, KGRPGL, XGRDGL ,YGRDGL, KGRBGL )
             !
             j = j + 1
             !
             xc(j) = xcl
             yc(j) = ycl
             !
          enddo
          !
       enddo
       !
       MXC    = itmp(1)
       MYC    = itmp(2)
       MCGRD  = itmp(3)
       NGRBND = itmp(4)
       MXF    = itmp(5)
       MYF    = itmp(6)
       !
       XCLMIN = rtmp(1)
       XCLMAX = rtmp(2)
       YCLMIN = rtmp(3)
       YCLMAX = rtmp(4)
       !
    endif
    !
    ! ... then carry out the interpolation for negative contribution to bulk dissipation
    !
    if ( optg /= 5 ) then
       !
       j = 0
       !
       do ixd = 1, mxd
          do iyd = 1, myd
             !
             j = j + 1
             !
             outsid = .false.
             !
             if ( xc(j) <= 0. ) then
                ixg = 1
                wi2 = 0.
                if ( xc(j) < -0.1 ) outsid = .true.
             else if ( xc(j) >= float(MXCGL-1) ) then
                ixg = MXCGL - 1
                wi2 = 1.
                if ( xc(j) > float(MXCGL)-0.9 ) outsid = .true.
             else
                ixg = int(xc(j))
                wi2 = xc(j) - real(ixg)
                ixg = ixg + 1
             endif
             !
             if ( yc(j) <= 0. ) then
                iyg = 1
                wj2 = 0.
                if ( yc(j) < -0.1 ) outsid = .true.
             else if ( yc(j) >= float(MYCGL-1) ) then
                iyg = MYCGL - 1
                wj2 = 1.
                if ( yc(j) > float(MYCGL)-0.9 ) outsid = .true.
             else
                iyg = int(yc(j))
                wj2 = yc(j) - real(iyg)
                iyg = iyg + 1
             endif
             !
             wi1 = 1.- wi2
             wj1 = 1.- wj2
             !
             if ( outsid ) then
                !
                dneg(iyd,ixd) = 0.
                !
             else
                !
                dneg(iyd,ixd) = wi1*wj1*disbk0(KGRPGL(ixg,iyg)) + wi1*wj2*disbk0(KGRPGL(ixg,iyg+1)) + wi2*wj1*disbk0(KGRPGL(ixg+1,iyg)) + wi2*wj2*disbk0(KGRPGL(ixg+1,iyg+1))
                !
             endif
             !
          enddo
       enddo
       !
    else
       !
       ! unstructured mesh
       !
       xp = xcugrd(vs(1))
       yp = ycugrd(vs(1))
       !
       do ixd = 1, mxd
          !
          x = xp - 0.5 * xpd(ixd)
          !
          do iyd = 1, myd
             !
             y = yp - 0.5 * ypd(iyd)
             !
             if ( x < XCGMIN .or. x > XCGMAX .or. y < YCGMIN .or. y > YCGMAX ) then
                !
                dneg(iyd,ixd) = 0.
                !
             else
                !
                call SwanInterpolatePoint( dneg(iyd,ixd), x, y, disbk0, 0. )
                !
             endif
             !
          enddo
          !
       enddo
       !
    endif
    !
    ! ... finally, compute the total bulk dissipation
    !
    disqc = 0.5 * ( dpos + dneg )
    !
    deallocate ( xc, yc, dpos, dneg )
    !
    end subroutine bdiss
    !
    subroutine varchk
    !
    use swcomm2
    use m_parall
    use SwanCompdata
    !
    integer, parameter :: maxit   = 2    ! maximum number of iterations
    real   , parameter :: larger  = 1.01 ! factor to make matrix element larger
    real   , parameter :: smaller = 0.99 ! factor to make matrix element smaller
    !
    integer            :: j              ! iteration counter
    real               :: rdev           ! indicates deviation in total dissipation in current grid point
    !
    j = 0
    !
 10 if ( varD > 0. .and. j < maxit ) then
       !
       j = j + 1
       !
       do ixd = 1, mxd
          do iyd = 1, myd
             !
             if ( disqc(iyd,ixd) > 0. ) then
                disqc(iyd,ixd) = smaller * disqc(iyd,ixd)
             elseif ( disqc(iyd,ixd) < 0. ) then
                disqc(iyd,ixd) = larger * disqc(iyd,ixd)
             endif
             !
          enddo
       enddo
       !
       varD = sum(disqc) * dkx * dky
       !
       goto 10
       !
    endif
    !
    if ( ITEST >= 50 ) then
       !
       if ( j > 0 ) then
          !
          ! check convergence
          !
          if ( optg /= 5 ) then
             write(PRTEST,'(a,i2,a,i5,a,i5,a)') ' ++ number of iterations to correct variance of bulk dissipation is ',j,' in grid point (',IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
          else
             write(PRTEST,'(a,i2,a,i7)') ' ++ number of iterations to correct variance of bulk dissipation is ',j,' in vertex k = ',vs(1)
          endif
          !
       elseif ( disbk < 0. .and. varW > 0. ) then
          !
          ! otherwise check consistency variance (see Eq. 25)
          !
          if ( varD /= varW * disbk ) then
             rdev = 100. * abs(varD - varW * disbk) / abs(varW * disbk)
             if ( rdev > 1. ) then
                if ( optg /= 5 ) then
                   write (PRTEST, '(a,f5.1,a,i5,a,i5,a)') 'integral of QC dissipation is not consistent - deviation=',rdev,' % in grid point (',IXCGRD(1)+MXF-1,',',IYCGRD(1)+MYF-1,')'
                else
                   write (PRTEST, '(a,f5.1,a,i7)') 'integral of QC dissipation is not consistent - deviation=',rdev,' % in vertex k = ',vs(1)
                endif
             endif
          endif
          !
       endif
    endif
    !
    if( varD < 0. .and. varW > 0. ) then
       !
       disqc = disqc * varW * disbk / varD
       !
    else
       !
       disqc = 0.
       !
    endif
    !
    end subroutine varchk
    !
end subroutine SWQCSURF
!
subroutine FILQCM ( imatra, idcmin, idcmax, isstop, memqcm, memqcb, plqcs, plwbrk, redc0, dissc0 )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!   41.90: Marcel Zijlema
!   41.91: Marcel Zijlema
!
!   Updates
!
!   41.90,     June 2021: New subroutine
!   41.91, February 2022: extension QCM with surf breaking
!
!   Purpose
!
!   fills IMATRA with the QC source terms for every gridpoint per sweep direction
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use swcomm4
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                         :: isstop ! maximum frequency that is propagated within a sweep
    !
    integer, dimension(MSC), intent(in)         :: idcmax ! maximum frequency-dependent counter in directional space
    integer, dimension(MSC), intent(in)         :: idcmin ! minimum frequency-dependent counter in directional space
    !
    real, dimension(MDC,MSC,MDISP), intent(out) :: dissc0 ! dissipation coefficient as explicit part (meant for output)
    real, dimension(MDC,MSC)      , intent(out) :: imatra ! coefficients of right hand side of action balance equation
    real, dimension(MDC,MSC,MCGRD), intent(in)  :: memqcb ! auxiliary array to store results of QC surf breaking in full spectral space
    real, dimension(MDC,MSC,MCGRD), intent(in)  :: memqcm ! auxiliary array to store results of QC scattering in full spectral space
    real, dimension(MDC,MSC,NPTST), intent(out) :: plqcs  ! array containing the QC scattering source term for test output
    real, dimension(MDC,MSC,NPTST), intent(out) :: plwbrk ! array containing the QC surf breaking source term for test output
    real, dimension(MDC,MSC,MREDS), intent(out) :: redc0  ! explicit part of redistribution in present vertex for output purposes
!
!   Local variables
!
    integer       :: id       ! loop counter over direction bins
    integer       :: iddum    ! counter in directional space
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: is       ! loop counter over frequency bins
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'FILQCM')
    !
    ! QC scattering
    !
    if ( IQCM > 0 ) then
       !
       do is = 1, isstop
          !
          do iddum = idcmin(is), idcmax(is)
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             ! store the results in the array IMATRA
             ! if TESTFL store results in array for isoline plot
             !
             imatra(id,is) = imatra(id,is) + memqcm(id,is,KCGRD(1))
             if ( TESTFL ) plqcs(id,is,IPTST) = memqcm(id,is,KCGRD(1))
             redc0(id,is,4) = redc0(id,is,4) + memqcm(id,is,KCGRD(1))
             !
          enddo
          !
       enddo
       !
       if ( TESTFL .and. ITEST > 50 ) then
          write (PRTEST,101) idcmin(1), idcmax(1), MSC, isstop
          if ( ITEST > 100 ) then
             do is = 1, isstop
                do iddum = idcmin(is), idcmax(is)
                   id = mod ( iddum - 1 + MDC , MDC ) + 1
                   write (PRTEST,102) is, id, memqcm(id,is,KCGRD(1))
                enddo
             enddo
          endif
       endif
       !
    endif
    !
    ! QC surf breaking
    !
    if ( ISURF > 0 .and. IGEN == 4 ) then
       !
       do is = 1, isstop
          !
          do iddum = idcmin(is), idcmax(is)
             id = mod ( iddum - 1 + MDC , MDC ) + 1
             !
             ! store the results in the array IMATRA
             ! if TESTFL store results in array for isoline plot
             !
             imatra(id,is) = imatra(id,is) + memqcb(id,is,KCGRD(1))
             if ( TESTFL ) plwbrk(id,is,IPTST) = memqcb(id,is,KCGRD(1))
             dissc0(id,is,2) = dissc0(id,is,2) + memqcb(id,is,KCGRD(1))
             !
          enddo
          !
       enddo
       !
       if ( TESTFL .and. ITEST > 50 ) then
          write (PRTEST,101) idcmin(1), idcmax(1), MSC, isstop
          if ( ITEST > 100 ) then
             do is = 1, isstop
                do iddum = idcmin(is), idcmax(is)
                   id = mod ( iddum - 1 + MDC , MDC ) + 1
                   write (PRTEST,103) is, id, memqcb(id,is,KCGRD(1))
                enddo
             enddo
          endif
       endif
       !
    endif
    !
 101 format(' FILQCM: ID_MIN ID_MAX MSC ISTOP :',4i6)
 102 format(' FILQCM: IS ID MEMQCM()          :',2i6,e12.4)
 103 format(' FILQCM: IS ID MEMQCB()          :',2i6,e12.4)
    !
end subroutine FILQCM
!
subroutine tukeywin ( a, n )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
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
!    41.90: Gal Akrish, Pieter Smit and Marcel Zijlema
!
!   Updates
!
!    41.90, June 2021: New subroutine
!
!   Purpose
!
!   Smooths function to be Fourier transformed
!
!   Method
!
!   The smoothing is based on a tapered cosine (Tukey) window
!
!   Modules used
!
    use ocpcomm4, only: ltrace
    use swcomm3 , only: pi
!
    implicit none
!
!   Argument variables
!
    integer     , intent(in)                    :: n ! size of array
    real(kind=8), dimension(n,n), intent(inout) :: a ! array to be smoothed
!
!   Parameter variables
!
    real, parameter :: w = 0.2  ! taper
!                                 (note: must be in between 0 and 0.5)
!
!   Local variables
!
    integer, save   :: ient = 0 ! number of entries in this subroutine
    integer         :: ix       ! loop counter in x-direction
    integer         :: iy       ! loop counter in y-direction
    !
    real            :: f        ! window fraction
    real            :: xw       ! window in x-direction
    real            :: yw       ! window in y-direction
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'tukeywin')
    !
    f = w * real(n-1)
    !
    do iy = 1, n
       !
       if ( real(iy - 1) < f ) then
          !
          yw = 1. + cos( pi * (real(iy - 1)/f - 1.) )
          !
       elseif ( real(n - iy) < f ) then
          !
          yw = 1. + cos( pi * (real(n - iy)/f - 1.) )
          !
       else
          !
          yw = 2.
          !
       endif
       !
       do ix = 1, n
          !
          if ( real(ix - 1) < f ) then
             !
             xw = 1. + cos( pi * (real(ix - 1)/f - 1.) )
             !
          elseif ( real(n - ix) < f ) then
             !
             xw = 1. + cos( pi * (real(n - ix)/f - 1.) )
             !
          else
             !
             xw = 2.
             !
          endif
          !
          a(ix,iy) = 0.25 * a(ix,iy) * xw * yw
          !
       enddo
       !
    enddo
    !
end subroutine tukeywin
!
end module SwanQCM
