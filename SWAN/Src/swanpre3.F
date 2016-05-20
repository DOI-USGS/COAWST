#include "swancpp.h"

      MODULE interp_swan_mod
#if defined WRF_COUPLING || defined NESTING
!
!=======================================================================
!  This module contains some original functions as well as             !
!  some routines modified from a ROMS subroutine interpolate.F and     !
!  contains several all purpuse generic routines:                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    compute_angle  computes swan grid angle.
!    swan_ref_init  Initialize swan refinement.                        !
!    shindices      Finds model grid cell for any datum.               !
!    stry_range     Binary search of model grid cell for any datum.    !
!    sinside        Closed polygon datum search.                       !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds

      implicit none

      CONTAINS

# ifdef WRF_COUPLING
      SUBROUTINE compute_angler (ng)
!
!=======================================================================
!                                                                      !
!  This routine computes the anlge of the swan grid.                   !
!                                                                      !
!=======================================================================
!
      USE M_GENARR
      USE M_PARALL
      USE SWCOMM2
      USE SWCOMM4
!
! Imported.
!
      integer, intent(in) :: ng
!
! Locals.
!
      integer :: i, j, count
      real, allocatable :: xlon(:,:), ylat(:,:)
      real, allocatable :: xlonb(:,:), ylatb(:,:)
      real :: cff1, cff2, cff3, angler, DEG2RAD
! 
! Compute SWAN grid angles for rotation of winds from WRF.
!
      ALLOCATE(M_GENARR_MOD(ng)%CosAngler_G(MXCGL*MYCGL))
      CosAngler=>M_GENARR_MOD(ng)%CosAngler_G
      ALLOCATE(M_GENARR_MOD(ng)%SinAngler_G(MXCGL*MYCGL))
      SinAngler=>M_GENARR_MOD(ng)%SinAngler_G
      ALLOCATE(xlon(MXCGL,MYCGL))
      ALLOCATE(ylat(MXCGL,MYCGL))
      ALLOCATE(xlonb(MXCGL+1,MYCGL))
      ALLOCATE(ylatb(MXCGL+1,MYCGL))
!
! xlon and xlat are the locations at the center points.
! xlonb and ylatb are offset by 1 and used to compute the angle
! at the center points.
!
      DO j=1,MYCGL
        DO i=1,MXCGL
          cff1=XOFFS_G(ng)
          cff2=PARALL_MOD(ng)%XGRDGL_G(i,j)
          xlon(i,j)=cff1+cff2
          cff1=YOFFS_G(ng)
          cff2=PARALL_MOD(ng)%YGRDGL_G(i,j)
          ylat(i,j)=cff1+cff2
        END DO
      END DO
!
      DO i=2,MXCGL
        DO j=1,MYCGL
          xlonb(i,j)=0.5*(xlon(i-1,j)+xlon(i,j))
          ylatb(i,j)=0.5*(ylat(i-1,j)+ylat(i,j))
        END DO
      END DO
      DO j=1,MYCGL
        xlonb(1,j)=xlonb(2,j)
        xlonb(MXCGL+1,j)=xlonb(MXCGL,j)
        ylatb(1,j)=ylatb(2,j)
        ylatb(MXCGL+1,j)=ylatb(MXCGL,j)
      END DO
!
! Now compute the angles.
!
      IF (KSPHER.eq.0) THEN ! Cartesian
        DEG2RAD=1.
      ELSE
        DEG2RAD=3.14159/180.
      END IF
      count=0
      DO j=1,MYCGL
        DO i=1,MXCGL
          count=count+1
          cff1=(xlonb(i+1,j)-xlonb(i,j))*DEG2RAD
          cff2=(ylatb(i+1,j)-ylatb(i,j))*DEG2RAD
          cff3=COS(ylat(i,j)*DEG2RAD)
          IF (KSPHER.eq.0) THEN ! Cartesian
            angler=ATAN2(cff2,cff1)
          ELSE
            angler=ATAN2(cff2,cff1*cff3)
          END IF
          CosAngler(count)=COS(angler)
          SinAngler(count)=SIN(angler)
        END DO
      END DO
!
      DEALLOCATE(xlon, ylat, xlonb, ylatb)
      RETURN
      END SUBROUTINE compute_angler
# endif

# ifdef NESTING
      SUBROUTINE swan_ref_init (ng, ngp)
!
!=======================================================================
!                                                                      !
!  This routine determines the bounding cells for each child grid,     !
!  determines the indices of the child cell in the parent grid,        !
!  and associates a parent tile for each bounding point.               !
!                                                                      !
!=======================================================================
!
      USE M_PARALL
      USE M_SREFINED
      USE SWCOMM2
      USE SWCOMM3
      USE M_GENARR

      integer, intent(in) :: ng, ngp

! Locals.
      integer :: i, j, count, IX, IY, MyError
      integer :: Xmin, Xmax, Ymin, Ymax
      INTEGER IARRL(5), IARRC(5,1:NPROC)
      real, allocatable :: Ipos(:,:), Jpos(:,:), angler(:,:)
      real :: IJspv
      logical :: rectangular
!
!  Create arrays of boundary points for current ng grid because
!  it is a child grid.
!
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDX_G(Numspec(ng)))
      BOUNDHINDX =>REFINED_MOD(ng)%BOUNDHINDX_G
      BOUNDHINDX=0.
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDY_G(Numspec(ng)))
      BOUNDHINDY =>REFINED_MOD(ng)%BOUNDHINDY_G
      BOUNDHINDY=0.
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDPF_G(Numspec(ng)))
      BOUNDHINDPF =>REFINED_MOD(ng)%BOUNDHINDPF_G
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDPR_G(Numspec(ng)))
      BOUNDHINDPR =>REFINED_MOD(ng)%BOUNDHINDPR_G
      IF (NPROC.gt.1) THEN
        BOUNDHINDPF=0
        BOUNDHINDPR=0
      ELSE
        BOUNDHINDPF=1
        BOUNDHINDPR=1
      END IF
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDIX_G(Numspec(ng)))
      BOUNDHINDIX =>REFINED_MOD(ng)%BOUNDHINDIX_G
      BOUNDHINDIX=0
      ALLOCATE (REFINED_MOD(ng)%BOUNDHINDIY_G(Numspec(ng)))
      BOUNDHINDIY =>REFINED_MOD(ng)%BOUNDHINDIY_G
      BOUNDHINDIY=0
!
!  Compute number of boundary array points to be passed.
!
      ac2size(ng)=MSC*MDC*Numspec(ng)
!
!  Determine the horizontal indices of each child
!  point in the parent grid.
!
      IJspv=9999.
      rectangular=.false.
      allocate(Ipos(MXCGL,MYCGL),Jpos(MXCGL,MYCGL))
      allocate(angler(MXCGL_G(ngp),MYCGL_G(ngp)))
      angler=0.
!
      CALL shindices(ng, 1, MXCGL_G(ngp), 1, MYCGL_G(ngp),              &
     &               1, MXCGL_G(ngp), 1, MYCGL_G(ngp),                  &
     &               angler,                                            &
     &               PARALL_MOD(ngp)%XGRDGL_G+XOFFS_G(ngp),             &
     &               PARALL_MOD(ngp)%YGRDGL_G+YOFFS_G(ngp),             &
     &               1, MXCGL, 1, MYCGL,                                &
     &               1, MXCGL, 1, MYCGL,                                &
     &               XGRDGL+XOFFS_G(ng), YGRDGL+YOFFS_G(ng), Ipos, Jpos,&
     &               IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
!
!  Fill the boundary arrays with the hindices.
!
      count=0
      DO i=1,MXCGL
        count=count+1
        IX=i
        IY=1
        BOUNDHINDX(count)=Ipos(IX,IY)
        BOUNDHINDY(count)=Jpos(IX,IY)
        BOUNDHINDIX(count)=IX
        BOUNDHINDIY(count)=IY
      END DO
      DO i=1,MXCGL
        count=count+1
        IX=i
        IY=MYCGL
        BOUNDHINDX(count)=Ipos(IX,IY)
        BOUNDHINDY(count)=Jpos(IX,IY)
        BOUNDHINDIX(count)=IX
        BOUNDHINDIY(count)=IY
      END DO
      DO i=2,MYCGL-1
        count=count+1
        IX=1
        IY=i
        BOUNDHINDX(count)=Ipos(IX,IY)
        BOUNDHINDY(count)=Jpos(IX,IY)
        BOUNDHINDIX(count)=IX
        BOUNDHINDIY(count)=IY
      END DO
      DO i=2,MYCGL-1
        count=count+1
        IX=MXCGL
        IY=i
        BOUNDHINDX(count)=Ipos(IX,IY)
        BOUNDHINDY(count)=Jpos(IX,IY)
        BOUNDHINDIX(count)=IX
        BOUNDHINDIY(count)=IY
      END DO
!
!  Determine which parent tile will give AC2 for each child point.
!
      IARRL(1) = MXF_G(ngp)
      IARRL(2) = MXL_G(ngp)
      IARRL(3) = MYF_G(ngp)
      IARRL(4) = MYL_G(ngp)
      IARRL(5) = MCGRD_G(ngp)
      CALL SWGATHER (IARRC, 5*NPROC, IARRL, 5, SWINT )
!
      CALL MPI_BCAST(IARRC,5*NPROC,SWINT,0,                             &
     &               WAV_COMM_WORLD,MyError)

!
!  BOUNDHINDPF = which processor is going to fill the AC2 for this point.
!  BOUNDHINDPR = which processor is able to read this data to fill its AC2 array.
!
      IF (NPROC.gt.1) THEN
        DO i=1,Numspec(ng)
          IF (MXCGL_G(ngp).gt.MYCGL_G(ngp)) THEN
            DO j=1,NPROC
              IF (BOUNDHINDPF(i).eq.0) THEN
                IF ((BOUNDHINDX(i).ge.IARRC(1,j)).and.                  &
     &              (BOUNDHINDX(i).le.IARRC(2,j)-3)) THEN
                  BOUNDHINDPF(i)=j
                ENDIF
              ENDIF
            END DO
          ELSE
            DO j=1,NPROC
              IF (BOUNDHINDPF(i).eq.0) THEN
                IF ((BOUNDHINDY(i).ge.IARRC(3,j)).and.                    &
     &              (BOUNDHINDY(i).le.IARRC(4,j)-3)) THEN
                  BOUNDHINDPF(i)=j
                ENDIF
              ENDIF
            END DO
          ENDIF
          IF (MXCGL.gt.MYCGL) THEN
            IF ((BOUNDHINDIX(i).ge.MXF).and.                               &
     &          (BOUNDHINDIX(i).le.MXL)) THEN
              BOUNDHINDPR(i)=INODE
            ENDIF
          ELSE
            IF ((BOUNDHINDIY(i).ge.MYF).and.                               &
     &          (BOUNDHINDIY(i).le.MYL)) THEN
              BOUNDHINDPR(i)=INODE
            ENDIF
          ENDIF
        END DO
      END IF

      deallocate(Ipos,Jpos)

      RETURN
      END SUBROUTINE swan_ref_init

      SUBROUTINE shindices (ng, LBi, UBi, LBj, UBj,                     &
     &                     Is, Ie, Js, Je,                              &
     &                     angler, Xgrd, Ygrd,                          &
     &                     LBm, UBm, LBn, UBn,                          &
     &                     Ms, Me, Ns, Ne,                              &
     &                     Xpos, Ypos, Ipos, Jpos,                      &
     &                     IJspv, rectangular, Xmin, Xmax, Ymin, Ymax)
!
!=======================================================================
!                                                                      !
!  Given any geographical locations Xpos and Ypos, this routine finds  !
!  the corresponding array cell indices (Ipos, Jpos) of gridded  data  !
!  Xgrd and Ygrd containing each requested location. This indices are  !
!  usually used for interpolation.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     LBi         I-dimension Lower bound of gridded data.             !
!     UBi         I-dimension Upper bound of gridded data.             !
!     LBj         J-dimension Lower bound of gridded data.             !
!     UBj         J-dimension Upper bound of gridded data.             !
!     Is          Starting gridded data I-index to search.             !
!     Ie          Ending   gridded data I-index to search.             !
!     Js          Starting gridded data J-index to search.             !
!     Je          Ending   gridded data J-index to search.             !
!     angler      Gridded data angle between X-axis and true EAST      !
!                   (radians).                                         !
!     Xgrd        Gridded data X-locations (usually, longitude).       !
!     Ygrd        Gridded data Y-locations (usually, latitude).        !
!     LBm         I-dimension Lower bound of requested locations.      !
!     UBm         I-dimension Upper bound of requested locations.      !
!     LBn         J-dimension Lower bound of requested locations.      !
!     UBn         J-dimension Upper bound of requested locations.      !
!     Ms          Starting requested locations I-index to search.      !
!     Me          Ending   requested locations I-index to search.      !
!     Ns          Starting requested locations J-index to search.      !
!     Ne          Ending   requested locations J-index to search.      !
!     Xpos        Requested X-locations to process (usually longitude).!
!     Ypos        Requested Y-locations to process (usually latitude). !
!     IJspv       Unbounded special value to assign.                   !
!     rectangular Logical switch indicating that gridded data has a    !
!                   plaid distribution.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-cell index containing locations in data. !
!     Jpos       Fractional J-cell index containing locations in data. !
!     Xmin       minimal cell index in x-dir containing data           !
!     Xmax       maximum cell index in x-dir containing data           !
!     Ymin       minimal cell index in y-dir containing data           !
!     Ymax       maximum cell index in y-dir containing data           !
!                                                                      !
!  Calls:    Try_Range                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
      USE SWCOMM3
      USE SWCOMM4
!
!  Imported variable declarations.
!
      logical, intent(in) :: rectangular
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj, Is, Ie, Js, Je
      integer, intent(in) :: LBm, UBm, LBn, UBn, Ms, Me, Ns, Ne

      real, intent(in) :: IJspv
!
      real, intent(in) :: angler(LBi:UBi,LBj:UBj)
      real, intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real, intent(in) :: Ygrd(LBi:UBi,LBj:UBj)

      real, intent(in) :: Xpos(LBm:UBm,LBn:UBn)
      real, intent(in) :: Ypos(LBm:UBm,LBn:UBn)

      real, intent(out) :: Ipos(LBm:UBm,LBn:UBn)
      real, intent(out) :: Jpos(LBm:UBm,LBn:UBn)

      integer, intent(out), optional :: Xmin, Xmax, Ymin, Ymax
!
!  Local variable declarations.
!
      logical :: found, foundi, foundj, foundxy

      integer :: Imax, Imin, Jmax, Jmin, i, i0, j, j0, mp, np

      real :: aa2, ang, bb2, diag2, dxx, dyy, phi
      real :: xfac, xpp, yfac, ypp
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
      Xmin=-9999
      Xmax=-9999
      Ymin=-9999
      Ymax=-9999
      foundxy=.TRUE.
      DO np=Ns,Ne
        DO mp=Ms,Me
          Ipos(mp,np)=IJspv
          Jpos(mp,np)=IJspv
!
!  The gridded data has a plaid distribution so the search is trivial.
!
          IF (rectangular) THEN
            foundi=.FALSE.
            I_LOOP : DO i=LBi,UBi-1
              IF ((Xgrd(i  ,1).le.Xpos(mp,np)).and.                     &
     &            (Xgrd(i+1,1).gt.Xpos(mp,np))) THEN
                Imin=i
                foundi=.TRUE.
                EXIT I_LOOP
              END IF
            END DO I_LOOP
            foundj=.FALSE.
            J_LOOP : DO j=LBj,UBj-1
              IF ((Ygrd(1,j  ).le.Ypos(mp,np)).and.                     &
     &            (Ygrd(1,j+1).gt.Ypos(mp,np))) THEN
                Jmin=j
                foundj=.TRUE.
                EXIT J_LOOP
              END IF
            END DO J_LOOP
            found=foundi.and.foundj
!
!  Check each position to find if it falls inside the whole domain.
!  Once it is stablished that it inside, find the exact cell to which
!  it belongs by successively dividing the domain by a half (binary
!  search).
!
          ELSE
            found=stry_range(ng, LBi, UBi, LBj, UBj,                    &
     &                      Xgrd, Ygrd,                                 &
     &                      Is, Ie, Js, Je,                             &
     &                      Xpos(mp,np), Ypos(mp,np))
            IF (found) THEN
              Imin=Is
              Imax=Ie
              Jmin=Js
              Jmax=Je
              DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
                IF ((Imax-Imin).gt.1) THEN
                  i0=(Imin+Imax)/2
                  found=stry_range(ng, LBi, UBi, LBj, UBj,              &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, i0, Jmin, Jmax,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Imax=i0
                  ELSE
                    Imin=i0
                  END IF
                END IF
                IF ((Jmax-Jmin).gt.1) THEN
                  j0=(Jmin+Jmax)/2
                  found=stry_range(ng, LBi, UBi, LBj, UBj,              &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, Imax, Jmin, j0,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Jmax=j0
                  ELSE
                    Jmin=j0
                  END IF
                END IF
              END DO
              found=(Is.le.Imin).and.(Imin.le.Ie).and.                  &
     &              (Is.le.Imax).and.(Imax.le.Ie).and.                  &
     &              (Js.le.Jmin).and.(Jmin.le.Je).and.                  &
     &              (Js.le.Jmax).and.(Jmax.le.Je)
            END IF
          END IF
!
!  Knowing the correct cell, calculate the exact indices, accounting
!  for a possibly rotated grid.  If spherical, convert all positions
!  to meters first.
!
          IF (found) THEN
            IF (KSPHER.gt.0) THEN
              yfac=REARTH*DEGRAD
              xfac=yfac*COS(Ypos(mp,np)*DEGRAD)
              xpp=(Xpos(mp,np)-Xgrd(Imin,Jmin))*xfac
              ypp=(Ypos(mp,np)-Ygrd(Imin,Jmin))*yfac
            ELSE
              xfac=1.0_m8
              yfac=1.0_m8
              xpp=Xpos(mp,np)-Xgrd(Imin,Jmin)
              ypp=Ypos(mp,np)-Ygrd(Imin,Jmin)
            END IF
!
!  Use Law of Cosines to get cell parallelogram "shear" angle.
!
            diag2=((Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+      &
     &            ((Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            aa2=((Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))*yfac)**2
            bb2=((Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            phi=ASIN((diag2-aa2-bb2)/(2.0_m8*SQRT(aa2*bb2)))
!
!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.
!
            ang=angler(Imin,Jmin)
            dxx=xpp*COS(ang)+ypp*SIN(ang)
            dyy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.
!
            dxx=dxx+dyy*TAN(phi)
            dyy=dyy/COS(phi)
!
!  Scale with cell side lengths to translate into cell indices.
!
            dxx=MIN(MAX(0.0,dxx/SQRT(aa2)),1.0)
            dyy=MIN(MAX(0.0,dyy/SQRT(bb2)),1.0)
            Ipos(mp,np)=REAL(Imin)+dxx
            Jpos(mp,np)=REAL(Jmin)+dyy
!
! Set min and max indices
!
            IF (foundxy) THEN
              Xmin=mp
              Ymin=np
              foundxy=.FALSE.
            ENDIF
            Xmax=MAX(Xmax, mp)
            Ymax=MAX(Ymax, np)
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE shindices

      LOGICAL FUNCTION stry_range (ng, LBi, UBi, LBj, UBj, Xgrd, Ygrd,  &
     &                            Imin, Imax, Jmin, Jmax, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this   !
!  function finds if the point (Xo,Yo)  is inside the box defined by   !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return   !
!  logical switch  try_range=.TRUE.  if (Xo,Yo) is inside, otherwise   !
!  it will return false.                                               !
!                                                                      !
!  Calls:   inside                                                     !
!                                                                      !
!=======================================================================
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax

#  ifdef ASSUMED_SHAPE
      real, intent(in) :: Xgrd(LBi:,LBj:)
      real, intent(in) :: Ygrd(LBi:,LBj:)
#  else
      real, intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real, intent(in) :: Ygrd(LBi:UBi,LBj:UBj)
#  endif

      real, intent(in) :: Xo, Yo
!
!  Local variable declarations.
!
      integer ::  Nb, i, j, shft, ic

      real, dimension(2*(Jmax-Jmin+Imax-Imin)+1) :: Xb, Yb
!
!-----------------------------------------------------------------------
!  Define closed polygon.
!-----------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      DO i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      END DO
      shft=1-Jmin+Imax-Imin
      DO j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      END DO
      shft=1+Jmax-Jmin+2*Imax-Imin
      DO i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      END DO
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      DO j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      END DO
!
!-----------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!-----------------------------------------------------------------------
!
      stry_range=sinside(Nb, Xb, Yb, Xo, Yo)
      RETURN
      END FUNCTION stry_range

      LOGICAL FUNCTION sinside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb

      real, intent(in) :: Xo, Yo

#  ifdef ASSUMED_SHAPE
      real, intent(inout) :: Xb(:), Yb(:)
#  else
      real, intent(inout) :: Xb(Nb+1), Yb(Nb+1)
#  endif
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex

      real :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_m8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0_m8).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0_m8).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0_m8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_m8)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        sinside=.FALSE.
      ELSE
        sinside=.TRUE.
      END IF
      RETURN
      END FUNCTION sinside
# endif
#endif
      END MODULE interp_swan_mod
