#include "w3macros.h"
!/ ------------------------------------------------------------------- /
MODULE W3PRO1MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           H. L. Tolman            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         05-Jun-2018 |
  !/                  +-----------------------------------+
  !/
  !/    04-Feb-2000 : Origination                         ( version 2.00 )
  !/    28-Mar-2001 : Partial time step bug fix (proper   ( version 2.10 )
  !/                  ingest of boundaries).
  !/    02-Apr-2001 : Sub-grid obstructions.              ( version 2.10 )
  !/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
  !/    20-Dec-2004 : Multiple grid version.              ( version 3.06 )
  !/    07-Sep-2005 : Improved XY boundary conditions.    ( version 3.08 )
  !/    10-Jan-2007 : Clean-up FACVX/Y compute.           ( version 3.10 )
  !/    05-Mar-2008 : Added NEC sxf90 compiler directives
  !/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
  !/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
  !/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
  !/                  (W. E. Rogers & T. J. Campbell, NRL)
  !/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
  !/                  specify index closure for a grid.   ( version 3.14 )
  !/                  (T. J. Campbell, NRL)
  !/    29-May-2014 : Adding OMPH switch.                 ( version 5.02 )
  !/    08-May-2014 : Implement tripolar grid for first order propagation
  !/                  scheme                              ( version 5.03 )
  !/                  (W. E. Rogers, NRL)
  !/    05-Jun-2018 : Add DEBUG                           ( version 6.04 )
  !/
  !/    Copyright 2009-2014 National Weather Service (NWS),
  !/       National Oceanic and Atmospheric Administration.  All rights
  !/       reserved.  WAVEWATCH III is a trademark of the NWS.
  !/       No unauthorized use without permission.
  !/
  !  1. Purpose :
  !
  !     Bundles routines for first order propagation scheme in single
  !     module.
  !
  !  2. Variables and types :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3MAP1    Subr. Public   Set up auxiliary maps.
  !      W3XYP1    Subr. Public   First order spatial propagation.
  !      W3KTP1    Subr. Public   First order spectral propagation.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      DSEC21    Func. W3TIMEMD Time difference.
  !      STRACE    Subr. W3SERVMD Subroutine tracing.
  !     ----------------------------------------------------------------
  !
  !  5. Remarks :
  !
  !  6. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/Tn    Enable test output.
  !
  !  7. Source code :
  !
  !/ ------------------------------------------------------------------- /
CONTAINS
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3MAP1 ( MAPSTA )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         06-Dec-2010 |
    !/                  +-----------------------------------+
    !/
    !/    19-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
    !/    14-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
    !/    20-Dec-2004 : Multiple grid version.              ( version 3.06 )
    !/    10-Jan-2007 : Clean-up FACVX/Y compute.           ( version 3.10 )
    !/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
    !/                  specify index closure for a grid.   ( version 3.14 )
    !/                  (T. J. Campbell, NRL)
    !/
    !  1. Purpose :
    !
    !     Generate 'map' arrays for the first order upstream scheme.
    !
    !  2. Method :
    !
    !     See section 3.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       MAPSTA  I.A.   I   Status map.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !     See module documentation.
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3WAVE    Subr. W3WAVEMD Wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     ------------------------------------------------------
    !      1.   Initialize arrays.
    !      2.   Fill arrays.
    !      3.   Invert arrays.
    !     ------------------------------------------------------
    !
    !  9. Switches :
    !
    !     !/S   Enable subroutine tracing.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE W3GDATMD, ONLY: NTH, NSPEC, NX, NY, ICLOSE,                 &
         ICLOSE_NONE, ICLOSE_SMPL, ICLOSE_TRPL
    USE W3ADATMD, ONLY: IS0, IS2, FACVX, FACVY
    USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: MAPSTA(NY*NX)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER                 :: IX, IY, IXY, ISP, IXNEXT
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3MAP1')
#endif
    !

    ! 1.  Initialize x-y arrays ------------------------------------------ *
    !
    FACVX = 0.
    FACVY = 0.
    !
    ! 2.  Fill x-y arrays ------------------------------------------------ *
    !
    !.....FACVY
    DO IX=1, NX
      DO IY=1, NY-1
        IXY    = IY +(IX-1)*NY
        IF ( MAPSTA( IXY ) .NE. 0 ) FACVY(IXY) = FACVY(IXY) + 1.
        !.........next point : j+1 : increment IXY by 1
        IF ( MAPSTA(IXY+1) .NE. 0 ) FACVY(IXY) = FACVY(IXY) + 1.
      END DO
    END DO
    !
    !.....FACVY for IY=NY
    IF ( ICLOSE.EQ.ICLOSE_TRPL ) THEN
      IY=NY
      DO IX=1, NX
        IXY    = IY +(IX-1)*NY
        IF ( MAPSTA( IXY ) .NE. 0 ) FACVY(IXY) = FACVY(IXY) + 1.
        !...........next point: j+1: tripole: j==>j+1==>j and i==>ni-i+1
        IXNEXT=NX-IX+1
        IXY    = IY +(IXNEXT-1)*NY
        IF ( MAPSTA( IXY ) .NE. 0 ) FACVY(IXY) = FACVY(IXY) + 1.
      END DO
      !BGR: Adding the following lines to compute FACVX over all
      !      IX for IY=NY (this allows along-seam propagation).
      !      Located here since already inside "TRPL" if-block.
      !{
      DO IX=1, NX-1
        IXY    = IY +(IX-1)*NY
        IF ( MAPSTA( IXY ) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
        IF ( MAPSTA(IXY+NY) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
      END DO
      !}
    END IF
    !
    !.....FACVX
    DO IX=1, NX-1
      DO IY=2, NY-1
        IXY    = IY +(IX-1)*NY
        IF ( MAPSTA( IXY  ) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
        !.........next point : i+1 : increment IXY by NY
        IF ( MAPSTA(IXY+NY) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
      END DO
    END DO
    !
    !.....FACVX for IX=NX
    IF ( ICLOSE.NE.ICLOSE_NONE ) THEN
      DO IY=2, NY-1
        IXY    = IY +(NX-1)*NY
        IF ( MAPSTA(IXY) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
        !...........next point : i+1 : increment IXY by NY
        !...........IXY+NY=IY+(IX-1)*NY+NY = IY+IX*NY = IY+NX*NY ==> wrap to IY
        IF ( MAPSTA(IY ) .NE. 0 ) FACVX(IXY) = FACVX(IXY) + 1.
      END DO
    END IF
    !
    ! 3.  Invert x-y arrays ---------------------------------------------- *
    !
    DO IXY=1, NX*NY
      IF ( FACVX(IXY) .NE. 0. ) FACVX(IXY) = 1. / FACVX(IXY)
      IF ( FACVY(IXY) .NE. 0. ) FACVY(IXY) = 1. / FACVY(IXY)
    END DO
    !
    ! 4.  Fill theta arrays ---------------------------------------------- *
    !
    DO ISP=1, NSPEC
      IS2  (ISP) = ISP + 1
      IS0  (ISP) = ISP - 1
    END DO
    !
    DO ISP=NTH, NSPEC, NTH
      IS2(ISP) = IS2(ISP) - NTH
    END DO
    !
    DO ISP=1, NSPEC, NTH
      IS0(ISP) = IS0(ISP) + NTH
    END DO
    !
    RETURN
    !/
    !/ End of W3MAP1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3MAP1
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3XYP1 ( ISP, DTG, MAPSTA, FIELD, VGX, VGY )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         29-May-2014 |
    !/                  +-----------------------------------+
    !/
    !/    07-Jul-1998 : Final FORTRAN 77                    ( version 1.18 )
    !/    14-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
    !/    28-Mar-2001 : Partial time step bug fix.          ( version 2.10 )
    !/    02-Apr-2001 : Sub-grid obstructions.              ( version 2.10 )
    !/    26-Dec-2002 : Moving grid version.                ( version 3.02 )
    !/    20-Dec-2004 : Multiple grid version.              ( version 3.06 )
    !/    07-Sep-2005 : Improved XY boundary conditions.    ( version 3.08 )
    !/    05-Mar-2008 : Added NEC sxf90 compiler directives
    !/                  (Chris Bunney, UK Met Office)       ( version 3.13 )
    !/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
    !/                  (W. E. Rogers & T. J. Campbell, NRL)
    !/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
    !/                  specify index closure for a grid.   ( version 3.14 )
    !/                  (T. J. Campbell, NRL)
    !/    29-May-2014 : Adding OMPH switch.                 ( version 5.02 )
    !/
    !  1. Purpose :
    !
    !     Propagation in physical space for a given spectral component.
    !
    !  2. Method :
    !
    !     First order scheme with flux formulation.
    !     Curvilinear grid implementation: Fluxes are computed in index space
    !       and then transformed back into physical space.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       ISP     Int.   I   Number of spectral bin (IK-1)*NTH+ITH
    !       DTG     Real   I   Total time step.
    !       MAPSTA  I.A.   I   Grid point status map.
    !       FIELD   R.A.  I/O  Wave action spectral densities on full
    !                          grid.
    !       VGX/Y   Real   I   Speed of grid.
    !     ----------------------------------------------------------------
    !
    !     Local variables.
    !     ----------------------------------------------------------------
    !       NTLOC   Int.  Number of local steps.
    !       DTLOC   Real  Local propagation time step.
    !       VCX     R.A.  Propagation velocities in index space.
    !       VCY     R.A.
    !       CXTOT   R.A.  Propagation velocities in physical space.
    !       CYTOT   R.A.
    !       VFLX    R.A.  Discrete fluxes between grid points in index space.
    !       VFLY    R.A.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !     See module documentation.
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3WAVE    Subr. W3WAVEMD Wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  7. Remarks :
    !
    !     - The local work arrays are initialized on the first entry to
    !       the routine.
    !     - Curvilinear grid implementation. Variables FACX, FACY, CCOS, CSIN,
    !       CCURX, CCURY are not needed and have been removed.  FACX is accounted
    !       for as approriate in this subroutine.  FACX is also accounted for in
    !       the case of .NOT.FLCX.  Since FACX is removed, there is now a check for
    !       .NOT.FLCX in this subroutine.  In CFL calcs dx and dy are omitted,
    !       since dx=dy=1 in index space.  Curvilinear grid derivatives
    !       (DPDY, DQDX, etc.) and metric (GSQRT) are brought in via W3GDATMD.
    !     - Standard VCB calculation for Y is:
    !               VCB       = FACVY(IXY) * ( VCY2D(IY,IX) + VCY2D(IY+1,IX) )
    !       This is to calculate the flux VCY(IY+0.5). For the tripole grid,
    !       we cannot do it this way, since the sign of VCY flips as we jump
    !       over the seam. If we were to do it this way, VCY(IY) and VCY(IY+1)
    !       are two numbers of similar magnitude and opposite sign, so the
    !       average of the two gives something close to zero, so energy does
    !       not leave via VCY(IY+0.5). One alternative is:
    !               VCB       = VCY2D(IY,IX)
    !       Another alternative is :
    !               VCB       = FACVY(IXY) * ( VCY2D(IY,IX) - VCY2D(IY+1,IX) )
    !       Both appear to give correct results for ww3_tp2.13. We use the
    !       second alternative.
    !
    !  8. Structure :
    !
    !     ---------------------------------------
    !       1.  Preparations
    !         a Set constants
    !         b Initialize arrays
    !       2.  Calculate local discrete fluxes
    !       3.  Calculate propagation fluxes
    !       4.  Propagate
    !       5.  Update boundary conditions
    !     ---------------------------------------
    !
    !  9. Switches :
    !
    !     !/S   Enable subroutine tracing.
    !
    !     !/OMPH  Hybrid OpenMP directives.
    !
    !     !/T   Enable general test output.
    !     !/T1  Test output local fluxes (V)FX-YL.
    !     !/T2  Test output propagation fluxes (V)FLX-Y.
    !     !/T3  Test output propagation.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS
    !
    USE W3TIMEMD, ONLY: DSEC21
    !
    USE W3GDATMD, ONLY: NK, NTH, SIG, ECOS, ESIN, NX, NY, NSEA,     &
         MAPSF, DTCFL, ICLOSE, CLATS, FLCX, FLCY,    &
         ICLOSE_NONE, ICLOSE_SMPL, ICLOSE_TRPL,      &
         FLAGLL, DPDX, DPDY, DQDX, DQDY, GSQRT
    USE W3WDATMD, ONLY: TIME
    USE W3ADATMD, ONLY: CG, CX, CY, ATRNX, ATRNY, FACVX, FACVY
    USE W3IDATMD, ONLY: FLCUR
    USE W3ODATMD, ONLY: NDST, FLBPI, NBI, TBPI0, TBPIN, ISBPI,      &
         BBPI0, BBPIN, NDSE, IAPROC, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: ISP, MAPSTA(NY*NX)
    REAL, INTENT(IN)        :: DTG, VGX, VGY
    REAL, INTENT(INOUT)     :: FIELD(1-NY:NY*(NX+2))
    !/
    !/ ------------------------------------------------------------------ /
    !/ Local parameters
    !/
    INTEGER                 :: IK, ITH, NTLOC, ITLOC, ISEA, IXY,    &
         IY0, IX, IY, JXN, JXP, JYN, JYP,     &
         IBI, NYMAX
#ifdef W3_T3
    INTEGER                 ::  IXF, IYF
#endif
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    REAL                    :: CG0, CGL, CGA, CC, CGN
    REAL                    :: DTLOC,DTRAD, VCB
    REAL                    :: RD1, RD2
    REAL                    :: CP, CQ
#ifdef W3_T3
    REAL                    :: AOLD
#endif
    !/
    !/ Automatic work arrays
    !/
    REAL                    :: CXTOT2D(NY,NX)
    REAL                    :: CYTOT2D(NY,NX)
    REAL                    :: FLD2D(NY+1,NX+1)
    REAL                    :: VCX2D(NY,NX+1)
    REAL                    :: VCY2D(NY+1,NX)
    REAL                    :: VFLX2D(1:NY,0:NX)
    REAL                    :: VFLY2D(NY,NX)

    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3XYP1')
#endif
    !
    ! 1.  Preparations --------------------------------------------------- *

    ! 1.a Set constants
    !
    ITH    = 1 + MOD(ISP-1,NTH)
    IK     = 1 + (ISP-1)/NTH
    !
    CG0    = 0.575 * GRAV / SIG(1)
    CGL    = 0.575 * GRAV / SIG(IK)
    !
    IF ( FLCUR ) THEN
      CGA    = SQRT(MAXVAL((CGL*ECOS(ITH)+CX(1:NSEA))**2          &
           +(CGL*ESIN(ITH)+CY(1:NSEA))**2))
      CC     = SQRT(MAXVAL(CX(1:NSEA)**2+CY(1:NSEA)**2))
#ifdef W3_MGP
      CGA    = SQRT(MAXVAL((CGL*ECOS(ITH)+CX(1:NSEA)-VGX)**2      &
           +(CGL*ESIN(ITH)+CY(1:NSEA)-VGY)**2))
      CC     = SQRT(MAXVAL((CX(1:NSEA)-VGX)**2+(CY(1:NSEA)-VGY)**2))
#endif
    ELSE
      CGA    = CGL
#ifdef W3_MGP
      CGA    = SQRT((CGL*ECOS(ITH)-VGX)**2+(CGL*ESIN(ITH)-VGY)**2)
#endif
      CC     = 0.
    END IF
    !
    CGN    = 0.9999 * MAX ( CGA, CC, 0.001*CG0 )
    !
    NTLOC  = 1 + INT(DTG/(DTCFL*CG0/CGN))
    DTLOC  = DTG / REAL(NTLOC)
    DTRAD  = DTLOC
    IF ( FLAGLL ) DTRAD=DTRAD/(DERA*RADIUS)

    !
#ifdef W3_T
    WRITE (NDST,9000) NTLOC
    WRITE (NDST,9001) ISP, ITH, IK
#endif
    !
    ! ====================== Loop partial ================================ *
    !
    DO ITLOC=1, NTLOC
      !
      ! 1.b Initialize arrays
      !
#ifdef W3_T1
      WRITE (NDST,9010) ITLOC
#endif
      !
      VCX2D = 0.
      VCY2D = 0.
      CXTOT2D  = 0.
      CYTOT2D  = 0.
      FLD2D  = 0.
      VFLX2D  = 0.
      VFLY2D  = 0.
      !
      ! 2.  Calculate field and velocities --------------------------------- *
      !
      !     FIELD = A / CG * CLATS
      !     VCX   = COS*CG / CLATS
      !     VCY   = SIN*CG
      !
#ifdef W3_T1
      WRITE (NDST,9020)
#endif
      !
#ifdef W3_OMPH
      !$OMP PARALLEL DO PRIVATE (ISEA, IXY, IX, IY)
#endif
      !
      DO ISEA=1, NSEA
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IXY    = MAPSF(ISEA,3)

        FLD2D(IY,IX) = FIELD(IXY) / CG(IK,ISEA) * CLATS(ISEA)

        CXTOT2D(IY,IX) = ECOS(ITH) * CG(IK,ISEA) / CLATS(ISEA)
        CYTOT2D(IY,IX) = ESIN(ITH) * CG(IK,ISEA)
#ifdef W3_MGP
        CXTOT2D(IY,IX) = CXTOT2D(IY,IX) - VGX/CLATS(ISEA)
        CYTOT2D(IY,IX) = CYTOT2D(IY,IX) - VGY
#endif

#ifdef W3_T1
        WRITE (NDST,9021) ISEA, IXY, FLD2D(IY,IX), &
             CXTOT2D(IY,IX), CYTOT2D(IY,IX)
#endif
      END DO
      !
#ifdef W3_OMPH
      !$OMP END PARALLEL DO
#endif
      !
      IF ( FLCUR ) THEN
        DO ISEA=1, NSEA
          IX     = MAPSF(ISEA,1)
          IY     = MAPSF(ISEA,2)

          CXTOT2D(IY,IX) = CXTOT2D(IY,IX) + CX(ISEA)/CLATS(ISEA)
          CYTOT2D(IY,IX) = CYTOT2D(IY,IX) + CY(ISEA)

        END DO
      END IF

      IF ( FLCX ) THEN
        DO ISEA=1, NSEA
          IX     = MAPSF(ISEA,1)
          IY     = MAPSF(ISEA,2)
          CP=CXTOT2D(IY,IX)*DPDX(IY,IX)+CYTOT2D(IY,IX)*DPDY(IY,IX)
          VCX2D(IY,IX) = CP*DTRAD
        END DO
      ELSE
        VCX2D=0.0
      ENDIF

      IF ( FLCY ) THEN
        DO ISEA=1, NSEA
          IX     = MAPSF(ISEA,1)
          IY     = MAPSF(ISEA,2)
          CQ=CXTOT2D(IY,IX)*DQDX(IY,IX)+CYTOT2D(IY,IX)*DQDY(IY,IX)
          VCY2D(IY,IX) = CQ*DTRAD
        END DO
      ELSE
        VCY2D=0.0
      ENDIF

      ! Transform FIELD to index space, i.e. straightened space
      ! Bugfix: This is now done *before* adding the ghost row, so that ghost
      !   row will be in index space (bug applied only to global, irregular
      !   grids, so it did not apply to any test case that existed w/v4.18)
      FLD2D(1:NY,1:NX)=FLD2D(1:NY,1:NX)*GSQRT(1:NY,1:NX)

      !
      ! Deal with longitude closure by duplicating one row *to the right*
      !   in FIELD/FLD2D, VCX
      IF ( ICLOSE.NE.ICLOSE_NONE ) THEN
#ifdef W3_T1
        WRITE (NDST,9024)
#endif
        DO IY=1, NY
          FLD2D(IY,NX+1)=FLD2D(IY,1)
          VCX2D(IY,NX+1)=VCX2D(IY,1)
#ifdef W3_T1
          WRITE (NDST,9025) IY, FLD2D(IY,NX+1), VCX2D(IY,NX+1)
#endif
        END DO
      END IF

      ! Deal with tripole closure by duplicating one row *at the top*
      !   in FIELD/FLD2D, VCY
      IF ( ICLOSE.EQ.ICLOSE_TRPL ) THEN
        DO IX=1,NX
          !...........next point: j+1: tripole: j==>j+1==>j and i==>ni-i+1
          FLD2D(NY+1,IX)=FLD2D(NY,NX-IX+1)
          VCY2D(NY+1,IX)=VCY2D(NY,NX-IX+1)
        END DO
      END IF

      !
      ! 3.  Calculate propagation fluxes ----------------------------------- *
      !
      NYMAX=NY-1
      IF ( ICLOSE.EQ.ICLOSE_TRPL ) NYMAX=NY
      !
#ifdef W3_OMPH
      !$OMP PARALLEL DO PRIVATE (IX, IY, IXY, VCB)
#endif
      !
      DO IX=1, NX
        DO IY=1, NYMAX
          IXY    = IY +(IX-1)*NY
          VCB       = FACVX(IXY) * ( VCX2D(IY,IX) + VCX2D(IY,IX+1) )
          VFLX2D(IY,IX) = MAX ( VCB , 0. ) * FLD2D(IY,IX)          &
               + MIN ( VCB , 0. ) * FLD2D(IY,IX+1)
        END DO
      END DO
      !
#ifdef W3_OMPH
      !$OMP END PARALLEL DO
#endif
      !
      ! Deal with longitude closure by duplicating one row *to the left*
      !    in VFLX. Note that a similar action is not take for tripole grid,
      !    since tripole seam is only: IY=NY communicating with other points
      !    at IY=NY, not a case of IY=NY communicating with IY=1
      IF ( ICLOSE.NE.ICLOSE_NONE ) THEN
#ifdef W3_T2
        WRITE (NDST,9032)
#endif
        DO IY=1, NY
          VFLX2D(IY,0) = VFLX2D(IY,NX)
#ifdef W3_T2
          WRITE (NDST,9033) IY, VFLX2D(IY,0)
#endif
        END DO
      END IF
      !
#ifdef W3_OMPH
      !$OMP PARALLEL DO PRIVATE (IX, IY, IXY, VCB)
#endif
      !
      DO IX=1, NX
        DO IY=1, NY-1
          IXY    = IY +(IX-1)*NY
          VCB       = FACVY(IXY) * ( VCY2D(IY,IX) + VCY2D(IY+1,IX) )
          VFLY2D(IY,IX) = MAX ( VCB , 0. ) * FLD2D(IY,IX)          &
               + MIN ( VCB , 0. ) * FLD2D(IY+1,IX)
        END DO
      END DO
      !
#ifdef W3_OMPH
      !$OMP END PARALLEL DO
#endif
      !

      ! For tripole grid, include IY=NY in calculation. VCB is handled
      !    differently. See notes in Section "7. Remarks" above.
      IF ( ICLOSE.EQ.ICLOSE_TRPL ) THEN
        IY=NY
        !
#ifdef W3_OMPH
        !$OMP PARALLEL DO PRIVATE (IXY, VCB, IX)
#endif
        !
        DO IX=1, NX
          IXY    = IY +(IX-1)*NY
          VCB       = FACVY(IXY) * ( VCY2D(IY,IX) - VCY2D(IY+1,IX) )
          VFLY2D(IY,IX) = MAX ( VCB , 0. ) * FLD2D(IY,IX)          &
               + MIN ( VCB , 0. ) * FLD2D(IY+1,IX)
        END DO
        !
#ifdef W3_OMPH
        !$OMP END PARALLEL DO
#endif
        !
      END IF

      ! 4.  Propagate ------------------------------------------------------ *
      !
#ifdef W3_T3
      WRITE (NDST,9040)
#endif
      !
#ifdef W3_OMPH
      !$OMP PARALLEL DO PRIVATE (ISEA, IXY, JXN, JXP, JYN, JYP, IX, IY)
#endif
      !
      DO ISEA=1, NSEA
        !
        IX    = MAPSF(ISEA,1)
        IY    = MAPSF(ISEA,2)
        IXY   = MAPSF(ISEA,3)

#ifdef W3_T3
        AOLD   = FLD2D(IY,IX) * CG(IK,ISEA) / CLATS(ISEA)
#endif
        !
        IF (MAPSTA(IXY).EQ.1) THEN
          !
          IF ( VFLX2D(IY,IX-1) .GT. 0. ) THEN
            JXN   = -1
          ELSE
            JXN   =  0
          END IF
          IF ( VFLX2D(IY,IX) .LT. 0. ) THEN
            JXP   =  1
          ELSE
            JXP   =  0
          END IF
          IF ( VFLY2D(IY-1,IX) .GT. 0. ) THEN
            JYN   = -1
          ELSE
            JYN   =  0
          END IF
          IF ( VFLY2D(IY,IX) .LT. 0. ) THEN
            JYP   =  1
          ELSE
            JYP   =  0
          END IF
          !
          FLD2D(IY,IX) =   FLD2D(IY,IX)                            &
               + ATRNX(IXY,JXN) * VFLX2D(IY,IX-1)   &
               - ATRNX(IXY,JXP) * VFLX2D(IY,IX)     &
               + ATRNY(IXY,JYN) * VFLY2D(IY-1,IX)   &
               - ATRNY(IXY,JYP) * VFLY2D(IY,IX)

#ifdef W3_T3
          WRITE (NDST,9041) ISEA, IXY, IXY-NY, IXY-1,         &
               VFLX2D(IY,IX-1), VFLX2D(IY,IX), VFLY2D(IY-1,IX), &
               VFLY2D(IY,IX) , CG(IK,ISEA)/CLATS(ISEA),AOLD,    &
               FLD2D(IY,IX)
#endif
          !
          !
#ifdef W3_T3
          WRITE (NDST,9042) ISEA, MAPSTA(IXY), AOLD,FLD2D(IY,IX)
#endif
          !
        END IF !  IF (MAPSTA(IXY).EQ.1) THEN

        FLD2D(IY,IX) = CG(IK,ISEA) / CLATS(ISEA) * FLD2D(IY,IX)

      END DO !  DO ISEA=1, NSEA
      !
#ifdef W3_OMPH
      !$OMP END PARALLEL DO
#endif
      !

      ! Transform FIELD back to physical space, i.e. may be curvilinear
      FLD2D(1:NY,1:NX)=FLD2D(1:NY,1:NX)/GSQRT(1:NY,1:NX)
      !
      ! 5.  Update boundary conditions ------------------------------------- *
      !
      IF ( FLBPI ) THEN
        RD1    = DSEC21 ( TBPI0, TIME ) - DTG *                   &
             REAL(NTLOC-ITLOC)/REAL(NTLOC)
        RD2    = DSEC21 ( TBPI0, TBPIN )
        IF ( RD2 .GT. 0.001 ) THEN
          RD2    = MIN(1.,MAX(0.,RD1/RD2))
          RD1    = 1. - RD2
        ELSE
          RD1    = 0.
          RD2    = 1.
        END IF
        DO IBI=1, NBI
          IX    = MAPSF(ISBPI(IBI),1)
          IY    = MAPSF(ISBPI(IBI),2)
          FLD2D(IY,IX) = RD1*BBPI0(ISP,IBI) + RD2*BBPIN(ISP,IBI)
        END DO
      END IF
      !
      ! 6.  Put back in 1d shape ------------------------------------------- *
      !
      DO ISEA=1, NSEA
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IXY    = MAPSF(ISEA,3)
        FIELD(IXY) = FLD2D(IY,IX)
      END DO
      !
      ! ... End of partial time step loop
      !
    END DO !   DO ITLOC=1, NTLOC
    !
    RETURN
    !
    ! Formats
    !
#ifdef W3_T
9000 FORMAT (' TEST W3XYP1 : NTLOC :',I4)
9001 FORMAT (' TEST W3XYP1 : ISP, ITH, IK :',I8,2I4)
#endif
    !
#ifdef W3_T1
9010 FORMAT (' TEST W3XYP1 : INIT. VFX-YL, ITLOC =',I3)
9020 FORMAT (' TEST W3XYP1 : ISEA, IXY, FIELD, VCX, VCY')
9021 FORMAT ('           ',2I8,3E12.4)
9024 FORMAT (' TEST W3XYP1 : GLOBAL CLOSURE: IY, FIELD, VCX ')
9025 FORMAT ('               ',I4,2E12.4)
#endif
    !
#ifdef W3_T2
9032 FORMAT (' TEST W3XYP1 : CLOSE. : IY, VFLX')
9033 FORMAT ('            ',I4,E12.4)
#endif
    !
#ifdef W3_T3
9040 FORMAT (' TEST W3XYP1 : PROPAGATION '/                      &
         '      ISEA, IXY(3), , FLX(2), FLY(2), FAC, A(2)')
9041 FORMAT (2X,4I5,1X,4E10.3,1X,E10.3,1X,2E10.3)
9042 FORMAT (2X,I5,'( MAP = ',I2,' )',56X,2E10.3)
#endif
    !/
    !/ End of W3XYP1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3XYP1
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3KTP1 ( ISEA, FACTH, FACK, CTHG0, CG, WN, DEPTH,    &
       DDDX, DDDY, CX, CY, DCXDX, DCXDY, DCYDX,    &
       DCYDY, DCDX, DCDY, VA )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         20-Dec-2004 |
    !/                  +-----------------------------------+
    !/
    !/    29-Aug-1997 : Final FORTRAN 77                    ( version 1.18 )
    !/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
    !/    20-Dec-2004 : Multiple grid version.              ( version 3.06 )
    !/
    !  1. Purpose :
    !
    !     Propagation in spectral space.
    !
    !  2. Method :
    !
    !     First order scheme.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       ISEA    Int.   I   Number of sea point.
    !       FACTH/K Real   I   Factor in propagation velocity.
    !       CTHG0   Real   I   Factor in great circle refracftion term.
    !       CG      R.A.   I   Local group velocities.
    !       WN      R.A.   I   Local wavenumbers.
    !       DEPTH   Real   I   Depth.
    !       DDDx    Real   I   Depth gradients.
    !       CX/Y    Real   I   Current components.
    !       DCxDx   Real   I   Current gradients.
    !       DCDX-Y  Real   I   Phase speed gradients.
    !       VA      R.A.  I/O  Spectrum.
    !     ----------------------------------------------------------------
    !
    !     Local variables.
    !     ----------------------------------------------------------------
    !       DSDD    R.A.  Partial derivative of sigma for depth.
    !       FRK, FRG, FKC
    !               R.A.  Partial velocity terms.
    !       DWNI    R.A.  Inverse band width.
    !       CTH-WN  R.A.  Propagation velocities of local fluxes.
    !       FLTH-WN R.A.  Propagation fluxes.
    !       AA      R.A.  Extracted spectrum
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !     See module documentation.
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3WAVE    Subr. W3WAVEMD Wave model routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !       None.
    !
    !  8. Structure :
    !
    !     -----------------------------------------------------------------
    !       1.  Preparations
    !         a Calculate DSDD
    !         b Extract spectrum
    !       2.  Refraction velocities
    !         a Filter level depth reffraction.
    !         b Depth refratcion velocity.
    !         c Current refraction velocity.
    !       3.  Wavenumber shift velocities
    !         a Prepare directional arrays
    !         b Calcuate velocity.
    !       4.  Refraction
    !         a Discrete fluxes.
    !         b Propagation fluxes.
    !         c Refraction.
    !       5.  Wavenumber shifts.
    !         a Discrete fluxes.
    !         b Propagation fluxes.
    !         c Refraction.
    !     -----------------------------------------------------------------
    !
    !  9. Switches :
    !
    !     C/S   Enable subroutine tracing.
    !     C/T   Enable general test output.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DSIP, ECOS, ESIN, ES2, &
         ESC, EC2, FACHFA, MAPWN, FLCTH, FLCK, CTMAX
    USE W3ADATMD, ONLY: IS0, IS2
    USE W3IDATMD, ONLY: FLCUR
    USE W3ODATMD, ONLY: NDST
#ifdef W3_DEBUG
    USE W3ODATMD, only : IAPROC
#endif
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !/
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: ISEA
    REAL, INTENT(IN)        :: FACTH, FACK, CTHG0, CG(0:NK+1),      &
         WN(0:NK+1), DEPTH, DDDX, DDDY,       &
         CX, CY, DCXDX, DCXDY, DCYDX, DCYDY
    REAL, INTENT(IN)        :: DCDX(0:NK+1), DCDY(0:NK+1)
    REAL, INTENT(INOUT)     :: VA(NSPEC)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
    INTEGER                 :: ITH, IK, ISP, ITH0
    REAL                    :: FDDMAX, FDG, DCYX, DCXXYY, DCXY,     &
         DCXX, DCXYYX, DCYY, FKD, FKD0, CTHB, &
         CWNB
    REAL                    :: VCTH(NSPEC), VCWN(1-NTH:NSPEC+NTH),  &
         VAA(1-NTH:NSPEC+NTH), VFLTH(NSPEC),  &
         VFLWN(1-NTH:NSPEC), DSDD(0:NK+1),    &
         FRK(NK), FRG(NK), FKC(NTH), DWNI(NK)
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
    CALL STRACE (IENT, 'W3KTP1')
#endif
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_T
    WRITE (NDST,9000) FLCTH, FLCK, FACTH, FACK, CTMAX
    WRITE (NDST,9001) ISEA, DEPTH, CX, CY,                       &
         DDDX, DDDY, DCXDX, DCXDY, DCYDX, DCYDY
#endif
    !
    ! 1.  Preparations --------------------------------------------------- *
    ! 1.a Array with partial derivative of sigma versus depth
    !
    DO IK=0, NK+1
      IF ( DEPTH*WN(IK) .LT. 5. ) THEN
        DSDD(IK) = MAX ( 0. ,                                     &
             CG(IK)*WN(IK)-0.5*SIG(IK) ) / DEPTH
      ELSE
        DSDD(IK) = 0.
      END IF
    END DO
    !
#ifdef W3_T
    WRITE (NDST,9010)
    DO IK=1, NK+1
      WRITE (NDST,9011) IK, TPI/SIG(IK), TPI/WN(IK),             &
           CG(IK), DSDD(IK)
    END DO
#endif
    !
    ! 1.b Extract spectrum
    !
    DO ISP=1, NSPEC
      VAA(ISP) = VA(ISP)
    END DO
    !
    ! 2.  Refraction velocities ------------------------------------------ *
    !
    IF ( FLCTH ) THEN
      !
      ! 2.a Set slope filter for depth refraction
      !
      FDDMAX = 0.
      FDG    = FACTH * CTHG0
      !
      DO ITH=1, NTH
        FDDMAX = MAX ( FDDMAX , ABS (                             &
             ESIN(ITH)*DDDX - ECOS(ITH)*DDDY ) )
      END DO
      !
      DO IK=1, NK
        FRK(IK) = FACTH * DSDD(IK) / WN(IK)
        FRK(IK) = FRK(IK) / MAX ( 1. , FRK(IK)*FDDMAX/CTMAX )
        FRG(IK) = FDG * CG(IK)
      END DO
      !
      ! 2.b Depth refraction and great-circle propagation
      !
      DO ISP=1, NSPEC
        VCTH(ISP) = FRG(MAPWN(ISP)) * ECOS(ISP)                   &
             + FRK(MAPWN(ISP)) * ( ESIN(ISP)*DDDX - ECOS(ISP)*DDDY )
      END DO
      !
#ifdef W3_REFRX
      ! 2.c @C/@x refraction and great-circle propagation
      VCTH = 0.
      FRK  = 0.
      FDDMAX = 0.
      !
      DO ISP=1, NSPEC
        FDDMAX = MAX ( FDDMAX , ABS (                      &
             ESIN(ISP)*DCDX(MAPWN(ISP)) - ECOS(ISP)*DCDY(MAPWN(ISP)) ) )
      END DO
      !
      DO IK=1, NK
        FRK(IK) = FACTH * CG(IK) * WN(IK) / SIG(IK)
        FRK(IK) = FRK(IK) / MAX ( 1. , FRK(IK)*FDDMAX/CTMAX )
        FRG(IK) = FDG * CG(IK)
      END DO
      DO ISP=1, NSPEC
        VCTH(ISP) = FRG(MAPWN(ISP)) * ECOS(ISP)            &
             + FRK(MAPWN(ISP)) * ( ESIN(ISP)*DCDX(MAPWN(ISP)) &
             - ECOS(ISP)*DCDY(MAPWN(ISP)) )
      END DO
#endif
      !
      ! 2.d Current refraction
      !
      IF ( FLCUR ) THEN
        !
        DCYX   = FACTH *   DCYDX
        DCXXYY = FACTH * ( DCXDX - DCYDY )
        DCXY   = FACTH *   DCXDY
        !
        DO ISP=1, NSPEC
          VCTH(ISP) = VCTH(ISP) + ES2(ISP)*DCYX                 &
               + ESC(ISP)*DCXXYY - EC2(ISP)*DCXY
        END DO
        !
      END IF
      !
    END IF
    !
    ! 3.  Wavenumber shift velocities ------------------------------------ *
    !
    IF ( FLCK ) THEN
      !
      DCXX   =  - FACK *   DCXDX
      DCXYYX =  - FACK * ( DCXDY + DCYDX )
      DCYY   =  - FACK *   DCYDY
      FKD    =    FACK * ( CX*DDDX + CY*DDDY )
      !
      DO ITH=1, NTH
        FKC(ITH) = EC2(ITH)*DCXX +                                &
             ESC(ITH)*DCXYYX + ES2(ITH)*DCYY
      END DO
      !
      ISP    = -NTH
      DO IK=0, NK+1
        FKD0   = FKD / CG(IK) * DSDD(IK)
        DO ITH=1, NTH
          ISP    = ISP + 1
          VCWN(ISP) = FKD0 + WN(IK)*FKC(ITH)
        END DO
      END DO
      !
      ITH0   = NSPEC - NTH
      DO ITH=1, NTH
        VAA(ITH+NSPEC) = FACHFA * VAA(ITH+ITH0)
        VAA(ITH- NTH ) = 0.
      END DO
      !
      DO IK=1, NK
        DWNI(IK) = CG(IK) / DSIP(IK)
      END DO
      !
    END IF
    !
    ! 4.  Refraction ----------------------------------------------------- *
    !
    IF ( FLCTH ) THEN
      !
      ! 4.a Boundary velocity and fluxes
      !
      DO ISP=1, NSPEC
        CTHB       = 0.5 * ( VCTH(ISP) + VCTH(IS2(ISP)) )
        VFLTH(ISP) = MAX ( CTHB , 0. ) * VAA(ISP)                 &
             + MIN ( CTHB , 0. ) * VAA(IS2(ISP))
      END DO
      !
      ! 4.b Propagation
      !
      DO ISP=1, NSPEC
        VA(ISP) = VA(ISP) + VFLTH(IS0(ISP)) - VFLTH(ISP )
      END DO
      !
    END IF
    !
    ! 5.  Wavenumber shifts ---------------------------------------------- *
    !
    IF ( FLCK ) THEN
      !
      ! 5.a Boundary velocity and fluxes
      !
      DO ISP=1-NTH, NSPEC
        CWNB       = 0.5 * ( VCWN(ISP) + VCWN(ISP+NTH) )
        VFLWN(ISP) = MAX ( CWNB , 0. ) * VAA(  ISP  )             &
             + MIN ( CWNB , 0. ) * VAA(ISP+NTH)
      END DO
      !
      ! 5.c Propagation
      !
      DO ISP=1, NSPEC
        VA(ISP) = VA(ISP) + DWNI(MAPWN(ISP)) *                    &
             ( VFLWN(ISP-NTH) - VFLWN(ISP) )
      END DO
      !
    END IF
    !
    RETURN
    !
    ! Formats
    !
#ifdef W3_T
9000 FORMAT (' TEST W3KTP1 : FLCTH-K, FACTH-K, CTMAX  :',         &
         2L2,2E10.3,F7.3)
9001 FORMAT (' TEST W3KTP1 : LOCAL DATA :',I7,F7.1,2F6.2,1X,      &
         6E10.3)
9010 FORMAT (' TEST W3KTP1 : IK, T, L, CG, DSDD : ')
9011 FORMAT ('              ',I3,F7.2,F7.1,F7.2,E11.3)
#endif
    !/
    !/ End of W3KTP1 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3KTP1
  !/
  !/ End of module W3PRO1MD -------------------------------------------- /
  !/
END MODULE W3PRO1MD
