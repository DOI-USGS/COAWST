!
! CRTM_Interpolation
!
! Module containing the generic polynomial interpolation
! routines used in the CRTM
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 01-Feb-2007
!                       paul.vandelst@noaa.gov
!
MODULE CRTM_Interpolation

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Type_Kinds, ONLY: fp
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Parameters
  PUBLIC :: ORDER
  PUBLIC :: NPTS
  ! Derived types and associated procedures
  PUBLIC :: LPoly_type
  PUBLIC :: Clear_LPoly
  PUBLIC :: LPoly_Init
  PUBLIC :: LPoly_Inspect
  ! Procedures
  PUBLIC :: Interp_1D
  PUBLIC :: Interp_2D
  PUBLIC :: Interp_3D
  PUBLIC :: Interp_4D
  PUBLIC :: Interp_1D_TL
  PUBLIC :: Interp_2D_TL
  PUBLIC :: Interp_3D_TL
  PUBLIC :: Interp_4D_TL
  PUBLIC :: Interp_1D_AD
  PUBLIC :: Interp_2D_AD
  PUBLIC :: Interp_3D_AD
  PUBLIC :: Interp_4D_AD
  PUBLIC :: Find_Index
  PUBLIC :: LPoly
  PUBLIC :: LPoly_TL
  PUBLIC :: LPoly_AD


  ! -------------------
  ! Procedure overloads
  ! -------------------
  INTERFACE Find_Index
    MODULE PROCEDURE Find_Regular_Index
    MODULE PROCEDURE Find_Random_Index
  END INTERFACE Find_Index


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_RCS_ID=&
  '$Id: CRTM_Interpolation.f90 60152 2015-08-13 19:19:13Z paul.vandelst@noaa.gov $'
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  INTEGER,  PARAMETER :: ORDER     = 2            ! Quadratic
  INTEGER,  PARAMETER :: NPOLY_PTS = ORDER+1      ! No. of points in each polynomial
  INTEGER,  PARAMETER :: NPTS      = NPOLY_PTS+1  ! No. of points total


  ! -----------------------
  ! Derived type definition
  ! -----------------------
  TYPE :: LPoly_type
!    PRIVATE
    INTEGER :: Order=ORDER
    INTEGER :: nPts =NPOLY_PTS
    ! Left and right side polynomials
    REAL(fp) :: lp_left(NPOLY_PTS)  = ZERO
    REAL(fp) :: lp_right(NPOLY_PTS) = ZERO
    ! Left and right side weighting factors
    REAL(fp) :: w_left  = ZERO
    REAL(fp) :: w_right = ZERO
    ! Polynomial numerator differences
    REAL(fp) :: dxi_left(NPOLY_PTS)  = ZERO
    REAL(fp) :: dxi_right(NPOLY_PTS) = ZERO
    ! Polynomial denominator differences
    REAL(fp) :: dx_left(NPOLY_PTS)  = ZERO
    REAL(fp) :: dx_right(NPOLY_PTS) = ZERO
  END TYPE LPoly_type


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!
! NAME:
!       Clear_LPoly
!
! PURPOSE:
!       Subroutine to reinitialise the LPoly_type structure
!
! CALLING SEQUENCE:
!       CALL Clear_LPoly(p)
!
! OUTPUT ARGUMENTS:
!       p:       Reinitialised Lagrangian polynomial structure
!                UNITS:      N/A
!                TYPE:       TYPE(LPoly_type)
!                DIMENSION:  Scalar
!                ATTRIBUTES: INTENT(IN OUT)
!
! SIDE EFFECTS:
!      If the structure contains data on input, it is replaced with the
!      reinitialisation value.
!
!--------------------------------------------------------------------------------
  SUBROUTINE Clear_LPoly(p)
    TYPE(LPoly_type), INTENT(IN OUT) :: p
    p%Order     = ORDER
    p%nPts      = NPOLY_PTS
    p%lp_left   = ZERO
    p%lp_right  = ZERO
    p%w_left    = ZERO
    p%w_right   = ZERO
    p%dxi_left  = ZERO
    p%dxi_right = ZERO
    p%dx_left   = ZERO
    p%dx_right  = ZERO
  END SUBROUTINE Clear_LPoly


  SUBROUTINE LPoly_Init(self)
    TYPE(LPoly_type), INTENT(OUT) :: self
    self%Order = ORDER
  END SUBROUTINE LPoly_Init


  SUBROUTINE LPoly_Inspect( self )
    TYPE(LPoly_type), INTENT(IN) :: self
    WRITE(*,'(1x,"LPoly OBJECT")')
    WRITE(*,'(3x,"Order : ",i0)') self%Order
    WRITE(*,'(3x,"nPts  : ",i0)') self%nPts
    WRITE(*,'(3x,"Left-side :")')
    WRITE(*,'(5x,"Weighting factor : ",es13.6)') self%w_left
    WRITE(*,'(5x,"Polynomial :")')
    WRITE(*,'(5(1x,es13.6,:))') self%lp_left
    WRITE(*,'(5x,"Numerator differences :")')
    WRITE(*,'(5(1x,es13.6,:))') self%dxi_left
    WRITE(*,'(5x,"Denominator differences :")')
    WRITE(*,'(5(1x,es13.6,:))') self%dx_left
    WRITE(*,'(3x,"Right-side :")')
    WRITE(*,'(5x,"Weighting factor : ",es13.6)') self%w_right
    WRITE(*,'(5x,"Polynomial :")')
    WRITE(*,'(5(1x,es13.6,:))') self%lp_right
    WRITE(*,'(5x,"Numerator differences :")')
    WRITE(*,'(5(1x,es13.6,:))') self%dxi_right
    WRITE(*,'(5x,"Denominator differences :")')
    WRITE(*,'(5(1x,es13.6,:))') self%dx_right
  END SUBROUTINE LPoly_Inspect


!--------------------------------------------------------------------------------
!
! NAME:
!       Interp_1D
!       Interp_2D
!       Interp_3D
!       Interp_4D
!
! PURPOSE:
!       Subroutines to perform interpolation:
!         o  1-D for z=f(u)
!         o  2-D for z=f(u,v)
!         o  3-D for z=f(u,v,w)
!         o  4-D for z=f(u,v,w,x)
!
! CALLING SEQUENCE:
!       CALL Interp_1D(z, ulp, z_int)
!       CALL Interp_2D(z, ulp, vlp, z_int)
!       CALL Interp_3D(z, ulp, vlp, wlp, z_int)
!       CALL Interp_4D(z, ulp, vlp, wlp, xlp, z_int)
!
! INPUT ARGUMENTS:
!       z:                    Data to interpolate
!                             UNITS:      Variable
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Rank-1 (N_PTS), or
!                                         Rank-2 (N_PTS,NPTS), or
!                                         Rank-3 (N_PTS,NPTS,N_PTS), or
!                                         Rank-4 (N_PTS,NPTS,N_PTS,NPTS), or
!                             ATTRIBUTES: INTENT(IN)
!
!       ulp, vlp, wlp, xlp:   Interpolating polynomial structures for the
!                             respective dimensions from previous calls to
!                             the LPoly() subroutine
!                             UNITS:      N/A
!                             TYPE:       LPoly_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       z_int:                Interpolation result
!                             UNITS:      Same as z
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!      Output z_int argument has INTENT(IN OUT) to prevent default reinitialisation
!      purely for computational speed.
!
!--------------------------------------------------------------------------------
  ! 1-D routine
  SUBROUTINE Interp_1D(z, ulp, &  ! Input
                       z_int   )  ! Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp
    REAL(fp),         INTENT(IN OUT) :: z_int  ! INTENT(IN OUT) to preclude reinitialisation
    ! Perform interpolation
    z_int = ( ulp%w_left  * ( ulp%lp_left(1) *z(1) + &
                              ulp%lp_left(2) *z(2) + &
                              ulp%lp_left(3) *z(3) ) ) + &
            ( ulp%w_right * ( ulp%lp_right(1)*z(2) + &
                              ulp%lp_right(2)*z(3) + &
                              ulp%lp_right(3)*z(4) ) )
  END SUBROUTINE Interp_1D

  ! 2-D routine
  SUBROUTINE Interp_2D(z, ulp, vlp, &  ! Input
                       z_int        )  ! Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp
    REAL(fp),         INTENT(IN OUT) :: z_int  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS)
    ! Interpolate z in u dimension for all v
    DO i = 1, NPTS
      CALL Interp_1D(z(:,i),ulp,a(i))
    END DO
    ! Interpolate z in w dimension
    CALL Interp_1D(a,vlp,z_int)
  END SUBROUTINE Interp_2D

  ! 3-D routine
  SUBROUTINE Interp_3D(z, ulp, vlp, wlp, &  ! Input
                       z_int             )  ! Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp, wlp
    REAL(fp),         INTENT(IN OUT) :: z_int  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS)
    ! Interpolate z in u,v dimension for all w
    DO i = 1, NPTS
      CALL Interp_2D(z(:,:,i),ulp,vlp,a(i))
    END DO
    ! Interpolate a in w dimension
    CALL Interp_1D(a,wlp,z_int)
  END SUBROUTINE Interp_3D

  ! 4-D routine
  SUBROUTINE Interp_4D(z, ulp, vlp, wlp, xlp, &  ! Input
                       z_int                  )  ! Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:,:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp, wlp, xlp
    REAL(fp),         INTENT(IN OUT) :: z_int  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS)
    ! Interpolate z in u,v,w dimension for all x
    DO i = 1, NPTS
      CALL Interp_3D(z(:,:,:,i),ulp,vlp,wlp,a(i))
    END DO
    ! Interpolate a in x dimension
    CALL Interp_1D(a,xlp,z_int)
  END SUBROUTINE Interp_4D


!--------------------------------------------------------------------------------
!
! NAME:
!       Interp_1D_TL
!       Interp_2D_TL
!       Interp_3D_TL
!       Interp_4D_TL
!
! PURPOSE:
!       Subroutines to perform tangent-linear interpolation:
!         o  1-D for z=f(u)
!         o  2-D for z=f(u,v)
!         o  3-D for z=f(u,v,w)
!         o  4-D for z=f(u,v,w,x)
!
! CALLING SEQUENCE:
!       CALL Interp_1D_TL(z, ulp, z_TL, ulp_TL, z_int_TL)
!       CALL Interp_2D_TL(z, ulp, vlp, z_TL, ulp_TL, vlp_TL, z_int_TL)
!       CALL Interp_3D_TL(z, ulp, vlp, wlp, z_TL, ulp_TL, vlp_TL, wlp_TL, z_int_TL)
!       CALL Interp_4D_TL(z, ulp, vlp, wlp, xlp, z_TL, ulp_TL, vlp_TL, wlp_TL, xlp_TL, z_int_TL)
!
! INPUT ARGUMENTS:
!       z:                       Data to interpolate
!                                UNITS:      Variable
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Rank-1 (N_PTS), or
!                                            Rank-2 (N_PTS,NPTS), or
!                                            Rank-3 (N_PTS,NPTS,N_PTS), or
!                                            Rank-4 (N_PTS,NPTS,N_PTS,NPTS), or
!                                ATTRIBUTES: INTENT(IN)
!
!       ulp, vlp, wlp, xlp:      Interpolating polynomial structures for the
!                                respective dimensions from previous calls to
!                                the LPoly() subroutine
!                                UNITS:      N/A
!                                TYPE:       LPoly_type
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN)
!
!       z_TL:                    Tangent-linear interpolation data.
!                                UNITS:      Same as z
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Same rank as z input
!                                ATTRIBUTES: INTENT(IN)
!
!       ulp_TL, vlp_TL,          Tangent-linear interpolating polynomial
!       wlp_TL, xlp_TL:          structures for the respective dimensions
!                                from previous calls to the LPoly_TL() subroutine.
!                                UNITS:      N/A
!                                TYPE:       LPoly_type
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       z_int_TL:                Tangent-linear interpolation result
!                                UNITS:      Same as z
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!      Output z_int_TL argument has INTENT(IN OUT) to prevent default
!      reinitialisation purely for computational speed.
!
!--------------------------------------------------------------------------------
  ! 1-D routine
  SUBROUTINE Interp_1D_TL( z   , ulp   , &  ! FWD Input
                           z_TL, ulp_TL, &  ! TL  Input
                           z_int_TL      )  ! TL  Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp
    REAL(fp),         INTENT(IN)     :: z_TL(:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp_TL
    REAL(fp),         INTENT(IN OUT) :: z_int_TL  ! INTENT(IN OUT) to preclude reinitialisation
    ! Perform TL interpolation
    z_int_TL = ( ulp%w_left    * ( ulp%lp_left(1)    * z_TL(1) + &
                                   ulp_TL%lp_left(1) * z(1)    + &
                                   ulp%lp_left(2)    * z_TL(2) + &
                                   ulp_TL%lp_left(2) * z(2)    + &
                                   ulp%lp_left(3)    * z_TL(3) + &
                                   ulp_TL%lp_left(3) * z(3)    ) ) + &
               ( ulp_TL%w_left * ( ulp%lp_left(1)    * z(1)    + &
                                   ulp%lp_left(2)    * z(2)    + &
                                   ulp%lp_left(3)    * z(3)    ) ) + &

               ( ulp%w_right    * ( ulp%lp_right(1)    * z_TL(2) + &
                                    ulp_TL%lp_right(1) * z(2)    + &
                                    ulp%lp_right(2)    * z_TL(3) + &
                                    ulp_TL%lp_right(2) * z(3)    + &
                                    ulp%lp_right(3)    * z_TL(4) + &
                                    ulp_TL%lp_right(3) * z(4)    ) ) + &
               ( ulp_TL%w_right * ( ulp%lp_right(1)    * z(2)    + &
                                    ulp%lp_right(2)    * z(3)    + &
                                    ulp%lp_right(3)    * z(4)    ) )

  END SUBROUTINE Interp_1D_TL

  ! 2-D routine
  SUBROUTINE Interp_2D_TL( z,    ulp   , vlp   , &  ! FWD Input
                           z_TL, ulp_TL, vlp_TL, &  ! TL  Input
                           z_int_TL              )  ! TL  Output
    REAL(fp),         INTENT(IN)     :: z(:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp
    REAL(fp),         INTENT(IN)     :: z_TL(:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp_TL, vlp_TL
    REAL(fp),         INTENT(IN OUT) :: z_int_TL  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_TL(NPTS)
    ! Interpolate z in v dimension for all w
    DO i = 1, NPTS
      CALL Interp_1D(z(:,i),ulp,a(i))
      CALL Interp_1D_TL(z(:,i),ulp,z_TL(:,i),ulp_TL,a_TL(i))
    END DO
    ! Interpolate z in w dimension
    CALL Interp_1D_TL(a,vlp,a_TL,vlp_TL,z_int_TL)
  END SUBROUTINE Interp_2D_TL

  ! 3-D routine
  SUBROUTINE Interp_3D_TL( z   , ulp   , vlp   , wlp   , &  ! FWD Input
                           z_TL, ulp_TL, vlp_TL, wlp_TL, &  ! TL  Input
                           z_int_TL                      )  ! TL  Output
    ! Arguments
    REAL(fp),         INTENT(IN)  :: z(:,:,:)
    TYPE(LPoly_type), INTENT(IN)  :: ulp, vlp, wlp
    REAL(fp),         INTENT(IN)  :: z_TL(:,:,:)
    TYPE(LPoly_type), INTENT(IN)  :: ulp_TL, vlp_TL, wlp_TL
    REAL(fp),         INTENT(IN OUT) :: z_int_TL  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_TL(NPTS)
    ! Interpolate z in u,v dimension for all w
    DO i = 1, NPTS
      CALL Interp_2D(z(:,:,i),ulp,vlp,a(i))
      CALL Interp_2D_TL(z(:,:,i),ulp,vlp,z_TL(:,:,i),ulp_TL,vlp_TL,a_TL(i))
    END DO
    ! Interpolate a in w dimension
    CALL Interp_1D_TL(a,wlp,a_TL,wlp_TL,z_int_TL)
  END SUBROUTINE Interp_3D_TL

  ! 4-D routine
  SUBROUTINE Interp_4D_TL( z   , ulp   , vlp   , wlp   , xlp   , &  ! FWD Input
                           z_TL, ulp_TL, vlp_TL, wlp_TL, xlp_TL, &  ! TL  Input
                           z_int_TL                              )  ! TL  Output
    ! Arguments
    REAL(fp),         INTENT(IN)  :: z(:,:,:,:)
    TYPE(LPoly_type), INTENT(IN)  :: ulp, vlp, wlp, xlp
    REAL(fp),         INTENT(IN)  :: z_TL(:,:,:,:)
    TYPE(LPoly_type), INTENT(IN)  :: ulp_TL, vlp_TL, wlp_TL, xlp_TL
    REAL(fp),         INTENT(IN OUT) :: z_int_TL  ! INTENT(IN OUT) to preclude reinitialisation
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_TL(NPTS)
    ! Interpolate z in u,v,w dimension for all x
    DO i = 1, NPTS
      CALL Interp_3D(z(:,:,:,i),ulp,vlp,wlp,a(i))
      CALL Interp_3D_TL(z(:,:,:,i),ulp,vlp,wlp,z_TL(:,:,:,i),ulp_TL,vlp_TL,wlp_TL,a_TL(i))
    END DO
    ! Interpolate a in x dimension
    CALL Interp_1D_TL(a,xlp,a_TL,xlp_TL,z_int_TL)
  END SUBROUTINE Interp_4D_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       Interp_1D_AD
!       Interp_2D_AD
!       Interp_3D_AD
!       Interp_4D_AD
!
! PURPOSE:
!       Subroutines to perform adjoint interpolation:
!         o  1-D for z=f(u)
!         o  2-D for z=f(u,v)
!         o  3-D for z=f(u,v,w)
!         o  4-D for z=f(u,v,w,x)
!
! CALLING SEQUENCE:
!       CALL Interp_1D_AD( z   , ulp   , &  ! FWD Input
!                          z_int_AD    , &  ! AD  Input
!                          z_AD, ulp_AD  )  ! AD  Output
!
!       CALL Interp_2D_AD( z   , ulp   , vlp   , &  ! FWD Input
!                          z_int_AD            , &  ! AD  Input
!                          z_AD, ulp_AD, vlp_AD  )  ! AD  Output
!
!       CALL Interp_3D_AD( z   , ulp   , vlp   , wlp   , &  ! FWD Input
!                          z_int_AD                    , &  ! AD  Input
!                          z_AD, ulp_AD, vlp_AD, wlp_AD  )  ! AD  Output
!
!       CALL Interp_4D_AD( z   , ulp   , vlp   , wlp   , xlp   , &  ! FWD Input
!                          z_int_AD                            , &  ! AD  Input
!                          z_AD, ulp_AD, vlp_AD, wlp_AD, xlp_AD  )  ! AD  Output
!
! INPUT ARGUMENTS:
!       z:                       Data to interpolate
!                                UNITS:      Variable
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Rank-1 (N_PTS), or
!                                            Rank-2 (N_PTS,NPTS), or
!                                            Rank-3 (N_PTS,NPTS,N_PTS), or
!                                            Rank-4 (N_PTS,NPTS,N_PTS,NPTS), or
!                                ATTRIBUTES: INTENT(IN)
!
!       ulp, vlp, wlp, xlp:      Interpolating polynomial structures for the
!                                respective dimensions from previous calls to
!                                the LPoly() subroutine
!                                UNITS:      N/A
!                                TYPE:       LPoly_type
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN)
!
!       z_int_AD:                Adjoint interpolate.
!                                UNITS:      N/A
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUT ARGUMENTS
!       z_AD:                    Adjoint interpolation data. Subsequently passed
!                                into the LPoly_AD() subroutine.
!                                UNITS:      N/A
!                                TYPE:       REAL(fp)
!                                DIMENSION:  Same rank as z input
!                                ATTRIBUTES: INTENT(IN OUT)
!
!       ulp_AD, vlp_AD,          Adjoint interpolating polynomial
!       wlp_AD, xlp_AD:          structures for respective dimensions.
!                                Subsequently passed into the LPoly_AD()
!                                subroutine.
!                                UNITS:      N/A
!                                TYPE:       LPoly_type
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------
  ! 1-D routine
  SUBROUTINE Interp_1D_AD( z   , ulp   , &  ! FWD Input
                           z_int_AD    , &  ! AD  Input
                           z_AD, ulp_AD  )  ! AD  Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp
    REAL(fp),         INTENT(IN OUT) :: z_int_AD
    REAL(fp),         INTENT(IN OUT) :: z_AD(:)
    TYPE(LPoly_type), INTENT(IN OUT) :: ulp_AD
    ! Local variables
    REAL(fp) :: wl_z_int_AD, wr_z_int_AD
    ! Perform adjoint interpolation
    ulp_AD%w_right = ulp_AD%w_right + &
                     ( z_int_AD * ( ( ulp%lp_right(1) * z(2) ) + &
                                    ( ulp%lp_right(2) * z(3) ) + &
                                    ( ulp%lp_right(3) * z(4) ) ) )

    ulp_AD%w_left = ulp_AD%w_left + &
                    ( z_int_AD * ( ( ulp%lp_left(1) * z(1) ) + &
                                   ( ulp%lp_left(2) * z(2) ) + &
                                   ( ulp%lp_left(3) * z(3) ) ) )

    wr_z_int_AD = ulp%w_right * z_int_AD
    ulp_AD%lp_right(1) = ulp_AD%lp_right(1) + ( wr_z_int_AD * z(2) )
    ulp_AD%lp_right(2) = ulp_AD%lp_right(2) + ( wr_z_int_AD * z(3) )
    ulp_AD%lp_right(3) = ulp_AD%lp_right(3) + ( wr_z_int_AD * z(4) )

    wl_z_int_AD = ulp%w_left * z_int_AD
    ulp_AD%lp_left(1) = ulp_AD%lp_left(1) + ( wl_z_int_AD * z(1) )
    ulp_AD%lp_left(2) = ulp_AD%lp_left(2) + ( wl_z_int_AD * z(2) )
    ulp_AD%lp_left(3) = ulp_AD%lp_left(3) + ( wl_z_int_AD * z(3) )

    z_AD(1) = z_AD(1) + ( wl_z_int_AD * ulp%lp_left(1) )

    z_AD(2) = z_AD(2) + ( wr_z_int_AD * ulp%lp_right(1) ) + &
                        ( wl_z_int_AD * ulp%lp_left(2)  )

    z_AD(3) = z_AD(3) + ( wr_z_int_AD * ulp%lp_right(2) ) + &
                        ( wl_z_int_AD * ulp%lp_left(3)  )

    z_AD(4) = z_AD(4) + ( wr_z_int_AD * ulp%lp_right(3) )

  END SUBROUTINE Interp_1D_AD

  ! 2-D routine
  SUBROUTINE Interp_2D_AD( z   , ulp   , vlp   , &  ! FWD Input
                           z_int_AD            , &  ! AD  Input
                           z_AD, ulp_AD, vlp_AD  )  ! AD  Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp
    REAL(fp),         INTENT(IN OUT) :: z_int_AD
    REAL(fp),         INTENT(IN OUT) :: z_AD(:,:)
    TYPE(LPoly_type), INTENT(IN OUT) :: ulp_AD, vlp_AD
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_AD(NPTS)
    ! Forward calculations
    ! Interpolate z in u dimension for all v
    DO i = 1, NPTS
      CALL Interp_1D(z(:,i),ulp,a(i))
    END DO
    ! Adjoint calculations
    ! Initialize local AD variables
    a_AD = ZERO
    ! Adjoint of z interpolation in v dimension
    CALL Interp_1D_AD(a,vlp,z_int_AD,a_AD,vlp_AD)
    ! Adjoint of z interpolation in u dimension for all v
    DO i = 1, NPTS
      CALL Interp_1D_AD(z(:,i),ulp,a_AD(i),z_AD(:,i),ulp_AD)
    END DO
  END SUBROUTINE Interp_2D_AD

  ! 3-D routine
  SUBROUTINE Interp_3D_AD( z   , ulp   , vlp   , wlp   , &  ! FWD Input
                           z_int_AD                    , &  ! AD  Input
                           z_AD, ulp_AD, vlp_AD, wlp_AD  )  ! AD  Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp, wlp
    REAL(fp),         INTENT(IN OUT) :: z_int_AD
    REAL(fp),         INTENT(IN OUT) :: z_AD(:,:,:)
    TYPE(LPoly_type), INTENT(IN OUT) :: ulp_AD, vlp_AD, wlp_AD
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_AD(NPTS)

    ! Forward calculations
    ! Interpolate z in u and v dimension for all w
    DO i = 1, NPTS
      CALL Interp_2D(z(:,:,i),ulp,vlp,a(i))
    END DO

    ! Adjoint calculations
    ! Initialize local AD variables
    a_AD = ZERO
    ! Adjoint of a interpolation in w dimension
    CALL Interp_1D_AD(a,wlp,z_int_AD,a_AD,wlp_AD)
    ! Adjoint of z interpolation in u and v dimension for all w
    DO i = 1, NPTS
      CALL Interp_2D_AD(z(:,:,i),ulp,vlp,a_AD(i),z_AD(:,:,i),ulp_AD,vlp_AD)
    END DO
  END SUBROUTINE Interp_3D_AD

  ! 4-D routine
  SUBROUTINE Interp_4D_AD( z   , ulp   , vlp   , wlp   , xlp   , &  ! FWD Input
                           z_int_AD                            , &  ! AD  Input
                           z_AD, ulp_AD, vlp_AD, wlp_AD, xlp_AD  )  ! AD  Output
    ! Arguments
    REAL(fp),         INTENT(IN)     :: z(:,:,:,:)
    TYPE(LPoly_type), INTENT(IN)     :: ulp, vlp, wlp, xlp
    REAL(fp),         INTENT(IN OUT) :: z_int_AD
    REAL(fp),         INTENT(IN OUT) :: z_AD(:,:,:,:)
    TYPE(LPoly_type), INTENT(IN OUT) :: ulp_AD, vlp_AD, wlp_AD, xlp_AD
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: a(NPTS), a_AD(NPTS)

    ! Forward calculations
    ! Interpolate z in u,v,w dimension for all x
    DO i = 1, NPTS
      CALL Interp_3D(z(:,:,:,i),ulp,vlp,wlp,a(i))
    END DO

    ! Adjoint calculations
    ! Initialize local AD variables
    a_AD = ZERO
    ! Adjoint of a interpolation in x dimension
    CALL Interp_1D_AD(a,xlp,z_int_AD,a_AD,xlp_AD)
    ! Adjoint of z interpolation in u,v,w dimension for all x
    DO i = 1, NPTS
      CALL Interp_3D_AD(z(:,:,:,i),ulp,vlp,wlp,a_AD(i),z_AD(:,:,:,i),ulp_AD,vlp_AD,wlp_AD)
    END DO
  END SUBROUTINE Interp_4D_AD


!--------------------------------------------------------------------------------
!
! NAME:
!       Find_Index
!
! PURPOSE:
!       Subroutines to search abscissa data for 4-pt interplation indices.
!
! CALLING SEQUENCE:
!       For regularly spaced x-data:
!         CALL Find_Index( x, dx,        &  ! Input
!                          x_int,        &  ! In/Output
!                          i1, i2,       &  ! Output
!                          out_of_bounds )  ! Output
!
!       For irregularly spaced x-data:
!         CALL Find_Index( x,            &  ! Input
!                          x_int,        &  ! In/Output
!                          i1, i2,       &  ! Output
!                          out_of_bounds )  ! Output
!
! INPUTS:
!       x:             Abscissa data.
!                      UNITS:      Variable
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Rank-1 (N)
!                      ATTRIBUTES: INTENT(IN)
!
!       dx:            Abscissa data spacing for the regularly spaced case.
!                      UNITS:      Same as x.
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! INPUTS/OUTPUTS
!       x_int:         On input : Abscissa value at which an interpolate
!                                 is desired.
!                      On output: Valid abscissa value at which an interpolate
!                                 will be computed.
!                      If x_int < x(1),          then x_int = x(1) upon exit.
!                         x_int > x(N),          then x_int = x(N) upon exit
!                         x(1) <= x_int <= x(N), then x_int is unchanged upon exit.
!                      Also see the out_of_bounds output argument.
!                      UNITS:      Same as x.
!                      TYPE:       REAL(fp)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUTS
!       i1, i2:        Begin and end indices in the input x-array to
!                      use for the 4-pt interpolation at the value x_int.
!                      Three cases are possible for a x array of length N:
!                        Normal    : x(i1) < x(i1+1) <= x_int <= x(i1+2) < x(i1+3)
!                        Left edge : x_int < x(2);   then i1=1,   i2=4
!                        Right edge: x_int > x(N-1); then i1=N-3, i2=N
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
!       out_of_bounds: Logical variable that identifies if the interpolate point,
!                      x_int, is within the bounds of the search array, x.
!                      If    x(1) <= x_int <= x(n), out_of_bounds == .FALSE.
!                      Else                         out_of_bounds == .TRUE.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!      Output arguments have INTENT(IN OUT) to prevent default
!      reinitialisation purely for computational speed.
!
!--------------------------------------------------------------------------------
  ! Find indices for regular spacing
  SUBROUTINE Find_Regular_Index(x, dx, x_int, i1, i2, out_of_bounds)
    REAL(fp), INTENT(IN)     :: x(:)
    REAL(fp), INTENT(IN)     :: dx
    REAL(fp), INTENT(IN OUT) :: x_int
    INTEGER , INTENT(IN OUT) :: i1, i2
    LOGICAL , INTENT(IN OUT) :: out_of_bounds
    INTEGER :: n
    n = SIZE(x)
    out_of_bounds = .FALSE.
    IF ( x_int < x(1) .OR. x_int > x(n) ) out_of_bounds = .TRUE.
    i1 = FLOOR((x_int-x(1))/dx)+1-(NPOLY_PTS/2)
    i1 = MIN(MAX(i1,1),n-NPOLY_PTS)
    i2 = i1 + NPOLY_PTS
    IF (out_of_bounds .AND. i1==1) THEN
      x_int = x(1)
    ELSE IF (out_of_bounds .AND. i2==n) THEN
      x_int = x(n)
    END IF
  END SUBROUTINE Find_Regular_Index

  ! Find indices for random spacing.
  ! Assumption is that x(1) <= xInt <= x(n)
  ! (despite the MIN/MAX test)
  SUBROUTINE Find_Random_Index(x, x_int, i1, i2, out_of_bounds)
    REAL(fp), INTENT(IN)     :: x(:)
    REAL(fp), INTENT(IN OUT) :: x_int
    INTEGER , INTENT(IN OUT) :: i1, i2
    LOGICAL , INTENT(IN OUT) :: out_of_bounds
    INTEGER :: k, n
    n = SIZE(x)
    out_of_bounds = .FALSE.
    IF ( x_int < x(1) .OR. x_int > x(n) ) out_of_bounds = .TRUE.
    DO k=1,n
      IF (x_int <= x(k) ) EXIT
    END DO
    i1 = MIN(MAX(1,k-1-(NPOLY_PTS/2)),n-NPOLY_PTS)
    i2 = i1 + NPOLY_PTS
    IF (out_of_bounds .AND. i1==1) THEN
      x_int = x(1)
    ELSE IF (out_of_bounds .AND. i2==n) THEN
      x_int = x(n)
    END IF
  END SUBROUTINE Find_Random_Index


!--------------------------------------------------------------------------------
!
! NAME:
!       LPoly
!
! PURPOSE:
!       Subroutines to compute the Lagrangian polynomial terms for interpolation.
!
! CALLING SEQUENCE:
!       CALL LPoly( x, x_int, &  ! Input
!                   p         )  ! Output
!
!
! INPUT ARGUMENTS:
!       x:         4-pt abscissa data.
!                  UNITS:      Variable
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Rank-1 (4)
!                  ATTRIBUTES: INTENT(IN)
!
!       x_int:     Abscissa value at which an interpolate is desired.
!                  Typically, x(2) <= xInt <= x(3), except for the data edges.
!                  UNITS:      Same as x.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       p:         Lagrangian polynomial structure for subsequent use in the
!                  various interpolation subroutines.
!                  UNITS:      N/A
!                  TYPE:       TYPE(LPoly_type)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!      Output p argument is INTENT(IN OUT) to prevent default reinitialisation
!      purely for computational speed.
!
!--------------------------------------------------------------------------------

  SUBROUTINE LPoly(x, x_int, p)
    REAL(fp),         INTENT(IN)     :: x(:)  ! Input
    REAL(fp),         INTENT(IN)     :: x_int ! Input
    TYPE(LPoly_type), INTENT(IN OUT) :: p     ! Output. INTENT(IN OUT) to preclude reinitialisation

    ! Compute the numerator differences
    CALL Compute_dxi(x(1:3),x_int,p%dxi_left)
    CALL Compute_dxi(x(2:4),x_int,p%dxi_right)

    ! Compute the denominator differences
    CALL Compute_dx(x(1:3),p%dx_left)
    CALL Compute_dx(x(2:4),p%dx_right)

    ! Compute the quadratic polynomials
    CALL Compute_QPoly(p%dxi_left , p%dx_left , p%lp_left)
    CALL Compute_QPoly(p%dxi_right, p%dx_right, p%lp_right)

    ! Polynomial weights
    IF ( x_int < x(2) ) THEN
      p%w_right = ZERO
      p%w_left  = ONE
    ELSE IF ( x_int > x(3) ) THEN
      p%w_right = ONE
      p%w_left  = ZERO
    ELSE
      p%w_right = p%dxi_left(2) / (-p%dx_left(3))
      p%w_left  = ONE - p%w_right
    END IF
  END SUBROUTINE LPoly


!--------------------------------------------------------------------------------
!
! NAME:
!       LPoly_TL
!
! PURPOSE:
!       Subroutines to compute the tangent-linear Lagrangian polynomial terms
!       for interpolation.
!
! CALLING SEQUENCE:
!       CALL LPoly_TL( x   , x_int   , &  ! FWD Input
!                      p             , &  ! FWD Input
!                      x_TL, x_int_TL, &  ! TL  Input
!                      p_TL            )  ! TL  Output
!
!
! INPUT ARGUMENTS:
!       x:         4-pt abscissa data.
!                  UNITS:      Variable
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Rank-1 (4)
!                  ATTRIBUTES: INTENT(IN)
!
!       x_int:     Abscissa value at which an interpolate is desired.
!                  Typically, x(2) <= xInt <= x(3), except for the data edges.
!                  UNITS:      Same as x.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
!       p:         Lagrangian polynomial structure from previous call to
!                  the LPoly() subroutine.
!                  UNITS:      N/A
!                  TYPE:       TYPE(LPoly_type)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
!       x_TL:      Tangent-linear 4-pt abscissa data.
!                  UNITS:      Same as x.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Rank-1 (4)
!                  ATTRIBUTES: INTENT(IN)
!
!       x_int_TL:  Tangent-linear abscissa value at which an interpolate
!                  is desired.
!                  UNITS:      Same as x.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       p_TL:      Tangent-linear Lagrangian polynomial structure for subsequent
!                  use in the various tangent-linear interpolation subroutines.
!                  UNITS:      N/A
!                  TYPE:       TYPE(LPoly_type)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!      Output p_TL argument is INTENT(IN OUT) to prevent default reinitialisation
!      purely for computational speed.
!
!--------------------------------------------------------------------------------

  SUBROUTINE LPoly_TL(x   , x_int   , p , &
                      x_TL, x_int_TL, p_TL)
    REAL(fp),         INTENT(IN)     :: x(:)      ! FWD Input
    REAL(fp),         INTENT(IN)     :: x_int     ! FWD Input
    TYPE(LPoly_type), INTENT(IN)     :: p         ! FWD Input
    REAL(fp),         INTENT(IN)     :: x_TL(:)   ! TL  Input
    REAL(fp),         INTENT(IN)     :: x_int_TL  ! TL  Input
    TYPE(LPoly_type), INTENT(IN OUT) :: p_TL      ! TL  Output. INTENT(IN OUT) to preclude reinitialisation

    ! Compute the tangent-linear numerator differences
    CALL Compute_dxi_TL(x_TL(1:3),x_int_TL,p_TL%dxi_left)
    CALL Compute_dxi_TL(x_TL(2:4),x_int_TL,p_TL%dxi_right)

    ! Compute the tangent-linear denominator differences
    CALL Compute_dx_TL(x_TL(1:3),p_TL%dx_left)
    CALL Compute_dx_TL(x_TL(2:4),p_TL%dx_right)

    ! Compute the tangent-linear quadratic polynomials
    CALL Compute_QPoly_TL(p%dxi_left   , p%dx_left    , p%lp_left,  &
                          p_TL%dxi_left, p_TL%dx_left , p_TL%lp_left)
    CALL Compute_QPoly_TL(p%dxi_right   , p%dx_right   , p%lp_right,  &
                          p_TL%dxi_right, p_TL%dx_right, p_TL%lp_right)

    ! Polynomial weights
    IF ( x_int < x(2) .OR. x_int > x(3) ) THEN
      p_TL%w_right = ZERO
      p_TL%w_left  = ZERO
    ELSE
      p_TL%w_right = -( p_TL%dxi_left(2) + (p%w_right*p_TL%dx_left(3)) ) / p%dx_left(3)
      p_TL%w_left  = -p_TL%w_right
    END IF
  END SUBROUTINE LPoly_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       LPoly_AD
!
! PURPOSE:
!       Subroutines to compute the Lagrangian polynomial adjoint terms
!       for interpolation.
!
! CALLING SEQUENCE:
!       CALL LPoly_AD( x   , x_int   , &  ! FWD Input
!                      p             , &  ! FWD Input
!                      p_AD          , &  ! AD  Input
!                      x_AD, x_int_AD  )  ! AD Output
!
! INPUT ARGUMENTS:
!       x:         4-pt abscissa data.
!                  UNITS:      Variable
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Rank-1 (4)
!                  ATTRIBUTES: INTENT(IN)
!
!       x_int:     Abscissa value at which an interpolate is desired.
!                  Typically, x(2) <= xInt <= x(3), except for the data edges.
!                  UNITS:      Same as x.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
!       p:         Lagrangian polynomial structure from previous call to
!                  the LPoly() subroutine.
!                  UNITS:      N/A
!                  TYPE:       TYPE(LPoly_type)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
!       p_AD:      Adjoint Lagrangian polynomial structure from previous calls
!                  to the various adjoint interpolation subroutines.
!                  *Note*: Components are modified (zeroed out) in this routine.
!                  UNITS:      N/A
!                  TYPE:       TYPE(LPoly_type)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUT ARGUMENTS:
!       x_AD:      Adjoint 4-pt abscissa data.
!                  UNITS:      Variable
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Rank-1 (4)
!                  ATTRIBUTES: INTENT(IN OUT)
!
!       x_int_AD:  Adjoint abscissa value at which an interpolate is desired.
!                  UNITS:      Same as x_AD.
!                  TYPE:       REAL(fp)
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN OUT)
!
! SIDE EFFECTS:
!      The adjoint input argument, p_AD, is modified in this subroutine. Its
!      components are reinitialised (zeroed out).
!
!--------------------------------------------------------------------------------

  SUBROUTINE LPoly_AD(x   , x_int   , p , &
                      p_AD, x_AD, x_int_AD)
    REAL(fp),         INTENT(IN)     :: x(:)      ! FWD Input
    REAL(fp),         INTENT(IN)     :: x_int     ! FWD Input
    TYPE(LPoly_type), INTENT(IN)     :: p         ! FWD Input
    TYPE(LPoly_type), INTENT(IN OUT) :: p_AD      ! AD  Input
    REAL(fp),         INTENT(IN OUT) :: x_AD(:)   ! AD  Output
    REAL(fp),         INTENT(IN OUT) :: x_int_AD  ! AD  Output

    ! Polynomial weights
    IF ( x_int < x(2) .OR. x_int > x(3) ) THEN
      p_AD%w_right = ZERO
      p_AD%w_left  = ZERO
    ELSE
      p_AD%w_right = p_AD%w_right - p_AD%w_left
      p_AD%w_left  = ZERO
      p_AD%dx_left(3)  = p_AD%dx_left(3)  - (p%w_right*p_AD%w_right/p%dx_left(3))
      p_AD%dxi_left(2) = p_AD%dxi_left(2) - (p_AD%w_right/p%dx_left(3))
      p_AD%w_right = ZERO
    END IF

    ! "Right" side quadratic
    CALL Compute_QPoly_AD(p%dxi_right   , p%dx_right,  &
                          p%lp_right    , &
                          p_AD%lp_right , &
                          p_AD%dxi_right, p_AD%dx_right)

    ! "Left" side quadratic
    CALL Compute_QPoly_AD(p%dxi_left   , p%dx_left,  &
                          p%lp_left    , &
                          p_AD%lp_left , &
                          p_AD%dxi_left, p_AD%dx_left)

    ! Compute the adjoint denominator differences
    CALL Compute_dx_AD(p_AD%dx_right,x_AD(2:4))
    CALL Compute_dx_AD(p_AD%dx_left ,x_AD(1:3))

    ! Compute the adjoint numerator differences
    CALL Compute_dxi_AD(p_AD%dxi_right,x_AD(2:4),x_int_AD)
    CALL Compute_dxi_AD(p_AD%dxi_left ,x_AD(1:3),x_int_AD)
  END SUBROUTINE LPoly_AD


!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! ------------------------------------------------
  ! Subroutines to compute the quadratic polynomials
  ! ------------------------------------------------
  ! Forward model
  SUBROUTINE Compute_QPoly(dxi, dx, lp)
    REAL(fp), INTENT(IN)     :: dxi(:)  ! Input
    REAL(fp), INTENT(IN)     :: dx(:)   ! Input
    REAL(fp), INTENT(IN OUT) :: lp(:)   ! Output. INTENT(IN OUT) to preclude reinitialisation
    lp(1) = dxi(2)*dxi(3) / ( dx(1)*dx(2))
    lp(2) = dxi(1)*dxi(3) / (-dx(1)*dx(3))
    lp(3) = dxi(1)*dxi(2) / ( dx(2)*dx(3))
  END SUBROUTINE Compute_QPoly

  ! Tangent-linear model
  SUBROUTINE Compute_QPoly_TL(dxi   , dx   , lp,  &
                              dxi_TL, dx_TL, lp_TL)
    REAL(fp), INTENT(IN)     :: dxi(:)     ! FWD Input
    REAL(fp), INTENT(IN)     :: dx(:)      ! FWD Input
    REAL(fp), INTENT(IN)     :: lp(:)      ! FWD Input
    REAL(fp), INTENT(IN)     :: dxi_TL(:)  ! TL  Input
    REAL(fp), INTENT(IN)     :: dx_TL(:)   ! TL  Input
    REAL(fp), INTENT(IN OUT) :: lp_TL(:)   ! TL  Output. INTENT(IN OUT) to preclude reinitialisation
    lp_TL(1) = ( (dxi(3)*dxi_TL(2)      ) + &
                 (dxi(2)*dxi_TL(3)      ) - &
                 (dx(2) *dx_TL(1) *lp(1)) - &
                 (dx(1) *dx_TL(2) *lp(1))   ) / (dx(1)*dx(2))
    lp_TL(2) = ( (dxi(3)*dxi_TL(1)      ) + &
                 (dxi(1)*dxi_TL(3)      ) + &
                 (dx(3) *dx_TL(1) *lp(2)) + &
                 (dx(1) *dx_TL(3) *lp(2))   ) / (-dx(1)*dx(3))
    lp_TL(3) = ( (dxi(2)*dxi_TL(1)      ) + &
                 (dxi(1)*dxi_TL(2)      ) - &
                 (dx(3) *dx_TL(2) *lp(3)) - &
                 (dx(2) *dx_TL(3) *lp(3))   ) / (dx(2)*dx(3))
  END SUBROUTINE Compute_QPoly_TL

  ! Adjoint model
  SUBROUTINE Compute_QPoly_AD(dxi   , dx, &
                              lp    , &
                              lp_AD , &
                              dxi_AD, dx_AD)
    REAL(fp), INTENT(IN)     :: dxi(:)     ! FWD Input
    REAL(fp), INTENT(IN)     :: dx(:)      ! FWD Input
    REAL(fp), INTENT(IN)     :: lp(:)      ! FWD Input
    REAL(fp), INTENT(IN OUT) :: lp_AD(:)   ! AD  Input
    REAL(fp), INTENT(IN OUT) :: dxi_AD(:)  ! AD  Output
    REAL(fp), INTENT(IN OUT) :: dx_AD(:)   ! AD  Output
    REAL(fp) :: d
    ! Adjoint of lp(3)
    d = lp_AD(3)/(dx(2)*dx(3))
    dxi_AD(1) = dxi_AD(1) + d*dxi(2)
    dxi_AD(2) = dxi_AD(2) + d*dxi(1)
    dx_AD(2)  = dx_AD(2)  - d*dx(3)*lp(3)
    dx_AD(3)  = dx_AD(3)  - d*dx(2)*lp(3)
    lp_AD(3) = ZERO
    ! Adjoint of lp(2)
    d = lp_AD(2)/(-dx(1)*dx(3))
    dxi_AD(1) = dxi_AD(1) + d*dxi(3)
    dxi_AD(3) = dxi_AD(3) + d*dxi(1)
    dx_AD(2)  = dx_AD(2)  + d*dx(3)*lp(2)
    dx_AD(3)  = dx_AD(3)  + d*dx(2)*lp(2)
    lp_AD(2) = ZERO
    ! Adjoint of lp(1)
    d = lp_AD(1)/(dx(1)*dx(2))
    dxi_AD(2) = dxi_AD(2) + d*dxi(3)
    dxi_AD(3) = dxi_AD(3) + d*dxi(2)
    dx_AD(1)  = dx_AD(1)  - d*dx(2)*lp(1)
    dx_AD(2)  = dx_AD(2)  - d*dx(1)*lp(1)
    lp_AD(1) = ZERO
  END SUBROUTINE Compute_QPoly_AD



  ! -------------------------------------------------------------
  ! Subroutines to compute the polynomial denominator differences
  ! -------------------------------------------------------------
  ! Forward model
  SUBROUTINE Compute_dx(x,dx)
    REAL(fp), INTENT(IN)     :: x(:)   ! Input
    REAL(fp), INTENT(IN OUT) :: dx(:)  ! Output. INTENT(IN OUT) to preclude reinitialisation
    dx(1) = x(1)-x(2)
    dx(2) = x(1)-x(3)
    dx(3) = x(2)-x(3)
  END SUBROUTINE Compute_dx

  ! Tangent-linear model
  SUBROUTINE Compute_dx_TL(x_TL,dx_TL)
    REAL(fp), INTENT(IN)     :: x_TL(:)   ! TL Input
    REAL(fp), INTENT(IN OUT) :: dx_TL(:)  ! TL Output. INTENT(IN OUT) to preclude reinitialisation
    dx_TL(1) = x_TL(1)-x_TL(2)
    dx_TL(2) = x_TL(1)-x_TL(3)
    dx_TL(3) = x_TL(2)-x_TL(3)
  END SUBROUTINE Compute_dx_TL

  ! Adjoint model
  SUBROUTINE Compute_dx_AD(dx_AD,x_AD)
    REAL(fp), INTENT(IN OUT) :: dx_AD(:)  ! AD Input
    REAL(fp), INTENT(IN OUT) :: x_AD(:)   ! AD Output
    x_AD(3) = x_AD(3) - dx_AD(3)
    x_AD(2) = x_AD(2) + dx_AD(3)
    dx_AD(3) = ZERO
    x_AD(3) = x_AD(3) - dx_AD(2)
    x_AD(1) = x_AD(1) + dx_AD(2)
    dx_AD(2) = ZERO
    x_AD(2) = x_AD(2) - dx_AD(1)
    x_AD(1) = x_AD(1) + dx_AD(1)
    dx_AD(1) = ZERO
  END SUBROUTINE Compute_dx_AD


  ! -----------------------------------------------------------
  ! Subroutines to compute the polynomial numerator differences
  ! -----------------------------------------------------------
  ! Forward model
  SUBROUTINE Compute_dxi(x,xi,dxi)
    REAL(fp), INTENT(IN)     :: x(:)    ! Input
    REAL(fp), INTENT(IN)     :: xi      ! Input
    REAL(fp), INTENT(IN OUT) :: dxi(:)  ! Output. INTENT(IN OUT) to preclude reinitialisation
    dxi(1) = xi - x(1)
    dxi(2) = xi - x(2)
    dxi(3) = xi - x(3)
  END SUBROUTINE Compute_dxi

  ! Tangent-linear model
  SUBROUTINE Compute_dxi_TL(x_TL,xi_TL,dxi_TL)
    REAL(fp), INTENT(IN)     :: x_TL(:)    ! TL Input
    REAL(fp), INTENT(IN)     :: xi_TL      ! TL Input
    REAL(fp), INTENT(IN OUT) :: dxi_TL(:)  ! TL Output. INTENT(IN OUT) to preclude reinitialisation
    dxi_TL(1) = xi_TL - x_TL(1)
    dxi_TL(2) = xi_TL - x_TL(2)
    dxi_TL(3) = xi_TL - x_TL(3)
  END SUBROUTINE Compute_dxi_TL

  ! Adjoint model
  SUBROUTINE Compute_dxi_AD(dxi_AD,x_AD,xi_AD)
    REAL(fp), INTENT(IN OUT) :: dxi_AD(:)  ! AD Input
    REAL(fp), INTENT(IN OUT) :: x_AD(:)    ! AD Output
    REAL(fp), INTENT(IN OUT) :: xi_AD      ! AD Output
    x_AD(1)   = x_AD(1) - dxi_AD(1)
    xi_AD     = xi_AD   + dxi_AD(1)
    dxi_AD(1) = ZERO
    x_AD(2)   = x_AD(2) - dxi_AD(2)
    xi_AD     = xi_AD   + dxi_AD(2)
    dxi_AD(2) = ZERO
    x_AD(3)   = x_AD(3) - dxi_AD(3)
    xi_AD     = xi_AD   + dxi_AD(3)
    dxi_AD(3) = ZERO
  END SUBROUTINE Compute_dxi_AD

END MODULE CRTM_Interpolation

