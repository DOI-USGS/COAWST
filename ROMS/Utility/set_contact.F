#include "cppdefs.h"
#ifdef NESTING
      SUBROUTINE set_contact (ng, model)
!
!svn $Id: set_contact.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads nested grids contact points for each contact     !
!  region from input NetCDF file and allocates/initializes all the     !
!  contact region structures declared in "mod_nesting".                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_nesting
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical, dimension((Ngrids-1)*2) :: Lcoincident
      logical, dimension((Ngrids-1)*2) :: Lcomposite
      logical, dimension((Ngrids-1)*2) :: Lmosaic
      logical, dimension((Ngrids-1)*2) :: Lrefinement

      logical, dimension((Ngrids-1)*2) :: Linterpolate
!
      integer :: cr, dg, ibry, ic, ig, ip,  m, rg, vindex
      integer :: my_Ncontact, my_Ngrids, my_nLweights, my_nQweights

      integer, dimension(Ngrids) :: my_Lm, my_Mm
      integer, dimension(Ngrids) :: refine_factor

      integer, dimension((Ngrids-1)*2) :: NpointsR
      integer, dimension((Ngrids-1)*2) :: NpointsU
      integer, dimension((Ngrids-1)*2) :: NpointsV
!
      real(r8), allocatable :: Lweight(:,:)
# ifdef QUADRATIC_WEIGHTS
      real(r8), allocatable :: Qweight(:,:)
# endif
      real(r8), allocatable :: Xrg(:), Yrg(:)
      real(r8), allocatable :: angle(:)
      real(r8), allocatable :: dmde(:), dndx(:)
      real(r8), allocatable :: f(:)
      real(r8), allocatable :: h(:)
      real(r8), allocatable :: mask(:)
      real(r8), allocatable :: pm(:), pn(:)
!
!-----------------------------------------------------------------------
!  Read in nesting grids contact point packed information into local
!  arrays.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, NGCname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Check if the number of nested grids is correct.
!
      CALL netcdf_get_dim (ng, model, NGCname, 'Ngrids',                &
     &                     my_Ngrids)
      IF (exit_flag.ne.NoError) RETURN

      IF (my_Ngrids.ne.Ngrids) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'inconsistent parameter, Ngrids = ',        &
     &                      Ngrids, my_Ngrids
  10      FORMAT (/,' READ_CONTACT - ', a, i4, 2x, i4,                  &
     &            /,16x,'in input file:'2x,a)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Check number of contact regions, Ncontact = (Ngrids-1)*2.
!
      CALL netcdf_get_dim (ng, model, NGCname, 'Ncontact',              &
     &                     my_Ncontact)
      IF (exit_flag.ne.NoError) RETURN

      IF (my_Ncontact.ne.(Ngrids-1)*2) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'inconsistent parameter, Ncontact = ',      &
     &                      (Ngrids-1)*2, my_Ncontact
        END IF
        exit_flag=5
        RETURN
      END IF
      Ncontact=(Ngrids-1)*2
!
!  Check number of contact points for bilinear interpolation weights.
!
      CALL netcdf_get_dim (ng, model, NGCname, 'nLweights',             &
     &                     my_nLweights)
      IF (exit_flag.ne.NoError) RETURN

      IF (my_nLweights.ne.4) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'inconsistent parameter, nLweights = ',     &
     &                      4, my_nLweights
        END IF
        exit_flag=5
        RETURN
      END IF

# ifdef QUADRATIC_WEIGHTS
!
!  Check number of contact points for quadratic interpolation weights.
!
      CALL netcdf_get_dim (ng, model, NGCname, 'nQweights',             &
     &                     my_nQweights)
      IF (exit_flag.ne.NoError) RETURN

      IF (my_nQweights.ne.9) THEN
        IF (Master) THEN
          WRITE (stdout,10) 'inconsistent parameter, nQweights = ',     &
     &                      9, my_nQweights
        END IF
        exit_flag=5
        RETURN
      END IF
# endif
!
!  Get number of contact points in input NetCDF file.  It may include
!  or not contact point over land when land/sea masking is activated.
!
      CALL netcdf_get_dim (ng, model, NGCname, 'datum',                 &
     &                     NCdatum)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in grid dimensions and if the grid order is correct.  The order
!  of the grids is important in nesting.
!
      CALL netcdf_get_ivar (ng, model, NGCname, 'Lm', my_Lm)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_get_ivar (ng, model, NGCname, 'Mm', my_Mm)
      IF (exit_flag.ne.NoError) RETURN

      DO ig=1,Ngrids
        IF (my_Lm(ig).ne.Lm(ig)) THEN
          IF (Master) THEN
            WRITE (stdout,10) 'inconsistent grid order, Lm = ',         &
     &                        Lm(ig), my_Lm(ig)
          END IF
          exit_flag=5
          RETURN
        END IF
        IF (my_Mm(ig).ne.Mm(ig)) THEN
          IF (Master) THEN
            WRITE (stdout,10) 'inconsistent grid order, Mm = ',         &
     &                        Mm(ig), my_Mm(ig)
          END IF
          exit_flag=5
          RETURN
        END IF
      END DO
!
!  Read in nesting types logical switches, 1:Ngrids.
!
      CALL netcdf_get_lvar (ng, model, NGCname, 'coincident',           &
     &                      Lcoincident)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_get_lvar (ng, model, NGCname, 'composite',            &
     &                      Lcomposite)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_get_lvar (ng, model, NGCname, 'mosaic',               &
     &                      Lmosaic)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_get_lvar (ng, model, NGCname, 'refinement',           &
     &                      Lrefinement)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in refinement factor from donor grid, 1:Ngrids.
!
      CALL netcdf_get_ivar (ng, model, NGCname, 'refine_factor',        &
     &                      refine_factor)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in vertical interpolation at contact points switch, 1:Ncontact.
!
      CALL netcdf_get_lvar (ng, model, NGCname, 'interpolate',          &
     &                      Linterpolate)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in data donor and data receiver grid number, 1:Ncontact.
!
      IF (.not.allocated(donor_grid)) THEN
        allocate ( donor_grid((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'donor_grid',           &
     &                      donor_grid)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(receiver_grid)) THEN
        allocate ( receiver_grid((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'receiver_grid',        &
     &                      receiver_grid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in donor grid (I,J) indices at PSI points used to extract
!  refined grid. Values have a fill value of -999 if not applicable.
!
      IF (.not.allocated(I_left)) THEN
        allocate ( I_left(Ngrids) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'I_left',               &
     &                      I_left)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(I_right)) THEN
        allocate ( I_right(Ngrids) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'I_right',              &
     &                      I_right)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(J_bottom)) THEN
        allocate ( J_bottom(Ngrids) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'J_bottom',             &
     &                      J_bottom)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(J_top)) THEN
        allocate ( J_top(Ngrids) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'J_top',                &
     &                      J_top)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in starting and ending contact points (RHO-, U-, V-) index in
!  packed data vectors for each contact region, 1:Ncontact.
!
      IF (.not.allocated(NstrR)) THEN
        allocate ( NstrR((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NstrR', NstrR)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(NendR)) THEN
        allocate ( NendR((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NendR', NendR)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(NstrU)) THEN
        allocate ( NstrU((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NstrU', NstrU)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(NendU)) THEN
        allocate ( NendU((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NendU', NendU)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(NstrV)) THEN
        allocate ( NstrV((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NstrV', NstrV)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (.not.allocated(NendV)) THEN
        allocate ( NendV((Ngrids-1)*2) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'NendV', NendV)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact region number for each contact point, 1:NCdatum.
!
      IF (.not.allocated(contact_region)) THEN
        allocate ( contact_region(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'contact_region',       &
     &                      contact_region)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in flag (1:NCdatum) to determine if the contact point is at
!  receiver grid interior (zero value) or on the boundary (1:western,
!  2:southern, 3:eastern, 4:northern).
!
      IF (.not.allocated(on_boundary)) THEN
        allocate ( on_boundary(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'on_boundary',          &
     &                      on_boundary)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in I-left and J-bottom indices of donor cell containing contact
!  point, 1:NCdatum.
!
      IF (.not.allocated(Idg_cp)) THEN
        allocate ( Idg_cp(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'Idg', Idg_cp)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(Jdg_cp)) THEN
        allocate ( Jdg_cp(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'Jdg', Jdg_cp)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in I- and J-indices of receiver grid contact point, 1:NCdatum.
!
      IF (.not.allocated(Irg_cp)) THEN
        allocate ( Irg_cp(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'Irg', Irg_cp)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(Jrg_cp)) THEN
        allocate ( Jrg_cp(NCdatum) )
      END IF
      CALL netcdf_get_ivar (ng, model, NGCname, 'Jrg', Jrg_cp)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in X- and Y-location of receiver grid contact point, 1:NCdatum.
!
      IF (.not.allocated(Xrg)) THEN
        allocate ( Xrg(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'Xrg', Xrg)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(Yrg)) THEN
        allocate ( Yrg(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'Yrg', Yrg)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point horizontal linear interpolation weights,
!  1:NCdatum.
!
      IF (.not.allocated(Lweight)) THEN
        allocate ( Lweight(my_nLweights,NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'Lweight', Lweight)
      IF (exit_flag.ne.NoError) RETURN

# ifdef QUADRATIC_WEIGHTS
!
!  Read in contact point horizontal quadratic interpolation weights,
!  1:NCdatum.
!
      IF (.not.allocated(Qweight)) THEN
        allocate ( Qweight(my_nQweights,NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'Qweight', Qweight)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Read in contact point bathymetry at RHO-points, 1:NCdatum.
!
      IF (.not.allocated(h)) THEN
        allocate ( h(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'h', h)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point Coriolis parameter at RHO-points, 1:NCdatum.
!
      IF (.not.allocated(f)) THEN
        allocate ( f(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'f', f)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point curvilinar coordinates in XI- and ETA-direction
!  at RHO-points, 1:NCdatum.
!
      IF (.not.allocated(pm)) THEN
        allocate ( pm(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'pm', pm)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(pn)) THEN
        allocate ( pn(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'pn', pn)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point inverse metric factor in XI- and ETA-direction,
!  1:NCdatum.
!
      IF (.not.allocated(dndx)) THEN
        allocate ( dndx(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'dndx', dndx)
      IF (exit_flag.ne.NoError) RETURN

      IF (.not.allocated(dmde)) THEN
        allocate ( dmde(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'dmde', dmde)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point angle between XI-axis and EAST, 1:NCdatum.
!
      IF (.not.allocated(angle)) THEN
        allocate ( angle(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'angle', angle)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in contact point land/sea mask, 1:NCdatum.
!
      IF (.not.allocated(mask)) THEN
        allocate ( mask(NCdatum) )
      END IF
      CALL netcdf_get_fvar (ng, model, NGCname, 'mask', mask)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Unpack contact point data into nesting structures (type T_NGC).
!-----------------------------------------------------------------------
!
!  Determine number of contact points in each contact region.
!
      IF (.not.allocated(NCpoints)) THEN
        allocate ( NCpoints(Ncontact) )
      END IF
      DO cr=1,Ncontact
        NpointsR(cr)=NendR(cr)-NstrR(cr)+1
        NpointsU(cr)=NendU(cr)-NstrU(cr)+1
        NpointsV(cr)=NendV(cr)-NstrV(cr)+1
        NCpoints(cr)=NpointsR(cr)+NpointsU(cr)+NpointsV(cr)
      END DO
!
!  Allocate grid connectivity (type T_NGC) structures.
!
      allocate ( Rcontact(Ncontact) )
      allocate ( Ucontact(Ncontact) )
      allocate ( Vcontact(Ncontact) )
!
!  Allocate arrays in grid connectivity structure.
!
      DO cr=1,Ncontact
        dg=donor_grid(cr)
        rg=receiver_grid(cr)

        allocate ( Rcontact(cr) % Irg(NpointsR(cr)) )
        allocate ( Ucontact(cr) % Irg(NpointsU(cr)) )
        allocate ( Vcontact(cr) % Irg(NpointsV(cr)) )

        allocate ( Rcontact(cr) % Jrg(NpointsR(cr)) )
        allocate ( Ucontact(cr) % Jrg(NpointsU(cr)) )
        allocate ( Vcontact(cr) % Jrg(NpointsV(cr)) )

        allocate ( Rcontact(cr) % Idg(NpointsR(cr)) )
        allocate ( Ucontact(cr) % Idg(NpointsU(cr)) )
        allocate ( Vcontact(cr) % Idg(NpointsV(cr)) )

        allocate ( Rcontact(cr) % Jdg(NpointsR(cr)) )
        allocate ( Ucontact(cr) % Jdg(NpointsU(cr)) )
        allocate ( Vcontact(cr) % Jdg(NpointsV(cr)) )

# ifdef SOLVE3D
        allocate ( Rcontact(cr) % Kdg(N(dg),NpointsR(cr)) )
        allocate ( Ucontact(cr) % Kdg(N(dg),NpointsU(cr)) )
        allocate ( Vcontact(cr) % Kdg(N(dg),NpointsV(cr)) )
# endif

        allocate ( Rcontact(cr) % Lweight(4,NpointsR(cr)) )
        allocate ( Ucontact(cr) % Lweight(4,NpointsU(cr)) )
        allocate ( Vcontact(cr) % Lweight(4,NpointsV(cr)) )
# ifdef WET_DRY
        allocate ( Rcontact(cr) % LweightUnmasked(4,NpointsR(cr)) )
        allocate ( Ucontact(cr) % LweightUnmasked(4,NpointsU(cr)) )
        allocate ( Vcontact(cr) % LweightUnmasked(4,NpointsV(cr)) )
# endif

# ifdef QUADRATIC_WEIGHTS
        allocate ( Rcontact(cr) % Qweight(9,NpointsR(cr)) )
        allocate ( Ucontact(cr) % Qweight(9,NpointsU(cr)) )
        allocate ( Vcontact(cr) % Qweight(9,NpointsV(cr)) )
#  ifdef WET_DRY
        allocate ( Rcontact(cr) % QweightUnmasked(9,NpointsR(cr)) )
        allocate ( Ucontact(cr) % QweightUnmasked(9,NpointsU(cr)) )
        allocate ( Vcontact(cr) % QweightUnmasked(9,NpointsV(cr)) )
#  endif
# endif

# ifdef SOLVE3D
        allocate ( Rcontact(cr) % Vweight(2,N(dg),NpointsR(cr)) )
        allocate ( Ucontact(cr) % Vweight(2,N(dg),NpointsU(cr)) )
        allocate ( Vcontact(cr) % Vweight(2,N(dg),NpointsV(cr)) )
# endif
      END DO
!
!  Initialize grid connectivity structure.
!
      DO cr=1,Ncontact
        dg=donor_grid(cr)
        rg=receiver_grid(cr)

        Rcontact(cr) % coincident = Lcoincident(rg)
        Ucontact(cr) % coincident = Lcoincident(rg)
        Vcontact(cr) % coincident = Lcoincident(rg)

        Rcontact(cr) % interpolate = Linterpolate(rg)
        Ucontact(cr) % interpolate = Linterpolate(rg)
        Vcontact(cr) % interpolate = Linterpolate(rg)

        Rcontact(cr) % donor_grid = dg
        Ucontact(cr) % donor_grid = dg
        Vcontact(cr) % donor_grid = dg

        Rcontact(cr) % receiver_grid = rg
        Ucontact(cr) % receiver_grid = rg
        Vcontact(cr) % receiver_grid = rg

        Rcontact(cr) % Npoints = NpointsR(cr)
        Ucontact(cr) % Npoints = NpointsU(cr)
        Vcontact(cr) % Npoints = NpointsV(cr)

        DO m=1,NpointsR(cr)
          ip=m+NstrR(cr)-1
          Rcontact(cr) % Irg(m) = Irg_cp(ip)
          Rcontact(cr) % Jrg(m) = Jrg_cp(ip)
          Rcontact(cr) % Idg(m) = Idg_cp(ip)
          Rcontact(cr) % Jdg(m) = Jdg_cp(ip)
# ifdef WET_DRY
          Rcontact(cr) % LweightUnmasked(1,m) = Lweight(1,ip)
          Rcontact(cr) % LweightUnmasked(2,m) = Lweight(2,ip)
          Rcontact(cr) % LweightUnmasked(3,m) = Lweight(3,ip)
          Rcontact(cr) % LweightUnmasked(4,m) = Lweight(4,ip)
# endif
          Rcontact(cr) % Lweight(1,m) = Lweight(1,ip)
          Rcontact(cr) % Lweight(2,m) = Lweight(2,ip)
          Rcontact(cr) % Lweight(3,m) = Lweight(3,ip)
          Rcontact(cr) % Lweight(4,m) = Lweight(4,ip)
# ifdef QUADRATIC_WEIGHTS
#  ifdef WET_DRY
          Rcontact(cr) % QweightUnmasked(1,m) = Qweight(1,ip)
          Rcontact(cr) % QweightUnmasked(2,m) = Qweight(2,ip)
          Rcontact(cr) % QweightUnmasked(3,m) = Qweight(3,ip)
          Rcontact(cr) % QweightUnmasked(4,m) = Qweight(4,ip)
          Rcontact(cr) % QweightUnmasked(5,m) = Qweight(5,ip)
          Rcontact(cr) % QweightUnmasked(6,m) = Qweight(6,ip)
          Rcontact(cr) % QweightUnmasked(7,m) = Qweight(7,ip)
          Rcontact(cr) % QweightUnmasked(8,m) = Qweight(8,ip)
          Rcontact(cr) % QweightUnmasked(9,m) = Qweight(9,ip)
#  endif
          Rcontact(cr) % Qweight(1,m) = Qweight(1,ip)
          Rcontact(cr) % Qweight(2,m) = Qweight(2,ip)
          Rcontact(cr) % Qweight(3,m) = Qweight(3,ip)
          Rcontact(cr) % Qweight(4,m) = Qweight(4,ip)
          Rcontact(cr) % Qweight(5,m) = Qweight(5,ip)
          Rcontact(cr) % Qweight(6,m) = Qweight(6,ip)
          Rcontact(cr) % Qweight(7,m) = Qweight(7,ip)
          Rcontact(cr) % Qweight(8,m) = Qweight(8,ip)
          Rcontact(cr) % Qweight(9,m) = Qweight(9,ip)
# endif
        END DO

        DO m=1,NpointsU(cr)
          ip=m+NstrU(cr)-1
          Ucontact(cr) % Irg(m) = Irg_cp(ip)
          Ucontact(cr) % Jrg(m) = Jrg_cp(ip)
          Ucontact(cr) % Idg(m) = Idg_cp(ip)
          Ucontact(cr) % Jdg(m) = Jdg_cp(ip)
# ifdef WET_DRY
          Ucontact(cr) % LweightUnmasked(1,m) = Lweight(1,ip)
          Ucontact(cr) % LweightUnmasked(2,m) = Lweight(2,ip)
          Ucontact(cr) % LweightUnmasked(3,m) = Lweight(3,ip)
          Ucontact(cr) % LweightUnmasked(4,m) = Lweight(4,ip)
# endif
          Ucontact(cr) % Lweight(1,m) = Lweight(1,ip)
          Ucontact(cr) % Lweight(2,m) = Lweight(2,ip)
          Ucontact(cr) % Lweight(3,m) = Lweight(3,ip)
          Ucontact(cr) % Lweight(4,m) = Lweight(4,ip)
# ifdef QUADRATIC_WEIGHTS
#  ifdef WET_DRY
          Ucontact(cr) % QweightUnmasked(1,m) = Qweight(1,ip)
          Ucontact(cr) % QweightUnmasked(2,m) = Qweight(2,ip)
          Ucontact(cr) % QweightUnmasked(3,m) = Qweight(3,ip)
          Ucontact(cr) % QweightUnmasked(4,m) = Qweight(4,ip)
          Ucontact(cr) % QweightUnmasked(5,m) = Qweight(5,ip)
          Ucontact(cr) % QweightUnmasked(6,m) = Qweight(6,ip)
          Ucontact(cr) % QweightUnmasked(7,m) = Qweight(7,ip)
          Ucontact(cr) % QweightUnmasked(8,m) = Qweight(8,ip)
          Ucontact(cr) % QweightUnmasked(9,m) = Qweight(9,ip)
#  endif
          Ucontact(cr) % Qweight(1,m) = Qweight(1,ip)
          Ucontact(cr) % Qweight(2,m) = Qweight(2,ip)
          Ucontact(cr) % Qweight(3,m) = Qweight(3,ip)
          Ucontact(cr) % Qweight(4,m) = Qweight(4,ip)
          Ucontact(cr) % Qweight(5,m) = Qweight(5,ip)
          Ucontact(cr) % Qweight(6,m) = Qweight(6,ip)
          Ucontact(cr) % Qweight(7,m) = Qweight(7,ip)
          Ucontact(cr) % Qweight(8,m) = Qweight(8,ip)
          Ucontact(cr) % Qweight(9,m) = Qweight(9,ip)
# endif
        END DO

        DO m=1,NpointsV(cr)
          ip=m+NstrV(cr)-1
          Vcontact(cr) % Irg(m) = Irg_cp(ip)
          Vcontact(cr) % Jrg(m) = Jrg_cp(ip)
          Vcontact(cr) % Idg(m) = Idg_cp(ip)
          Vcontact(cr) % Jdg(m) = Jdg_cp(ip)
# ifdef WET_DRY
          Vcontact(cr) % LweightUnmasked(1,m) = Lweight(1,ip)
          Vcontact(cr) % LweightUnmasked(2,m) = Lweight(2,ip)
          Vcontact(cr) % LweightUnmasked(3,m) = Lweight(3,ip)
          Vcontact(cr) % LweightUnmasked(4,m) = Lweight(4,ip)
# endif
          Vcontact(cr) % Lweight(1,m) = Lweight(1,ip)
          Vcontact(cr) % Lweight(2,m) = Lweight(2,ip)
          Vcontact(cr) % Lweight(3,m) = Lweight(3,ip)
          Vcontact(cr) % Lweight(4,m) = Lweight(4,ip)
# ifdef QUADRATIC_WEIGHTS
#  ifdef WET_DRY
          Vcontact(cr) % QweightUnmasked(1,m) = Qweight(1,ip)
          Vcontact(cr) % QweightUnmasked(2,m) = Qweight(2,ip)
          Vcontact(cr) % QweightUnmasked(3,m) = Qweight(3,ip)
          Vcontact(cr) % QweightUnmasked(4,m) = Qweight(4,ip)
          Vcontact(cr) % QweightUnmasked(5,m) = Qweight(5,ip)
          Vcontact(cr) % QweightUnmasked(6,m) = Qweight(6,ip)
          Vcontact(cr) % QweightUnmasked(7,m) = Qweight(7,ip)
          Vcontact(cr) % QweightUnmasked(8,m) = Qweight(8,ip)
          Vcontact(cr) % QweightUnmasked(9,m) = Qweight(9,ip)
#  endif
          Vcontact(cr) % Qweight(1,m) = Qweight(1,ip)
          Vcontact(cr) % Qweight(2,m) = Qweight(2,ip)
          Vcontact(cr) % Qweight(3,m) = Qweight(3,ip)
          Vcontact(cr) % Qweight(4,m) = Qweight(4,ip)
          Vcontact(cr) % Qweight(5,m) = Qweight(5,ip)
          Vcontact(cr) % Qweight(6,m) = Qweight(6,ip)
          Vcontact(cr) % Qweight(7,m) = Qweight(7,ip)
          Vcontact(cr) % Qweight(8,m) = Qweight(8,ip)
          Vcontact(cr) % Qweight(9,m) = Qweight(9,ip)
# endif
        END DO

      END DO
!
!-----------------------------------------------------------------------
!  Load contact points grid metrics.
!-----------------------------------------------------------------------
!
!  Allocate contact region metrics (type T_NGM) structure.
!
      allocate ( CONTACT_METRIC(Ncontact) )
!
!  Allocate arrays in contact region metrics structure.
!
      DO cr=1,Ncontact
        allocate ( CONTACT_METRIC(cr) % angler(NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % dndx  (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % dmde  (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % f     (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % h     (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % rmask (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % umask (NpointsU(cr)) )
        allocate ( CONTACT_METRIC(cr) % vmask (NpointsV(cr)) )
        allocate ( CONTACT_METRIC(cr) % pm    (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % pn    (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % Xr    (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % Yr    (NpointsR(cr)) )
        allocate ( CONTACT_METRIC(cr) % Xu    (NpointsU(cr)) )
        allocate ( CONTACT_METRIC(cr) % Yu    (NpointsU(cr)) )
        allocate ( CONTACT_METRIC(cr) % Xv    (NpointsV(cr)) )
        allocate ( CONTACT_METRIC(cr) % Yv    (NpointsV(cr)) )
      END DO
!
!  Initialize contact region metrics structure.
!
      DO cr=1,Ncontact
        DO m=1,NpointsR(cr)
          ip=m+NstrR(cr)-1
          CONTACT_METRIC(cr) % angler(m) = angle(ip)
          CONTACT_METRIC(cr) % dndx  (m) = dndx (ip)
          CONTACT_METRIC(cr) % dmde  (m) = dmde (ip)
          CONTACT_METRIC(cr) % f     (m) = f    (ip)
          CONTACT_METRIC(cr) % h     (m) = h    (ip)
          CONTACT_METRIC(cr) % rmask (m) = mask (ip)
          CONTACT_METRIC(cr) % pm    (m) = pm   (ip)
          CONTACT_METRIC(cr) % pn    (m) = pn   (ip)
          CONTACT_METRIC(cr) % Xr    (m) = Xrg  (ip)
          CONTACT_METRIC(cr) % Yr    (m) = Yrg  (ip)
        END DO

        DO m=1,NpointsU(cr)
          ip=m+NstrU(cr)-1
          CONTACT_METRIC(cr) % umask(m) = mask(ip)
          CONTACT_METRIC(cr) % Xu   (m) = Xrg (ip)
          CONTACT_METRIC(cr) % Yu   (m) = Yrg (ip)
        END DO

        DO m=1,NpointsV(cr)
          ip=m+NstrV(cr)-1
          CONTACT_METRIC(cr) % vmask(m) = mask(ip)
          CONTACT_METRIC(cr) % Xv   (m) = Xrg (ip)
          CONTACT_METRIC(cr) % Yv   (m) = Yrg (ip)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initialize various parameters.
!-----------------------------------------------------------------------
!
!  Nested grids conectivity switches.  They are used in "get_bounds.F"
!  to determenine if the global state array need to be extended by few
!  extra points at each boundary to accomodate the contact regions and
!  contact points.  These are VERY important switches in ROMS nesting.
!
      IF (.not.allocated(ContactRegion)) THEN
        allocate ( ContactRegion(4,Ngrids) )
        ContactRegion = .FALSE.
      END IF
!
      DO m=1,NCdatum
        cr=contact_region(m)
        rg=receiver_grid(cr)
        ibry=on_boundary(m)
        IF ((ibry.eq.iwest ).or.(ibry.eq.ieast ).or.                    &
     &      (ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          IF (.not.ContactRegion(ibry,rg)) THEN
            ContactRegion(ibry,rg)=.TRUE.
          END IF
        END IF
      END DO
!
!  Set composite and refinement grids switches.
!
      DO cr=1,Ncontact
        rg=receiver_grid(cr)
        DO ibry=1,4
          IF (Lcomposite(cr).and.ContactRegion(ibry,rg)) THEN
            CompositeGrid(ibry,rg)=.TRUE.
          END IF
        END DO
        RefinedGrid(rg)=Lrefinement(cr)
      END DO
!
!  Refinement grid scale factor.
!
      DO ig=1,Ngrids
        RefineScale(ig)=Refine_factor(ig)
      END DO
!
!  Coarser donor grid number to finer grid contact points.  Their
!  values are computed in "metrics.F" when the maximum grid spacing
!  is computed. The maximum grid spacing is used to distinguish
!  contact areas for more than two layers of nesting.
!
!$OMP PARALLEL
      IF (.not.allocated(CoarserDonor)) THEN
        allocate ( CoarserDonor(Ngrids) )
        CoarserDonor = 0
      END IF
!$OMP END PARALLEL
!
!  Finer donor grid to coarser receiver grid for two-way feedback.
!  Their values are computed in "metrics.F" when the maximum grid
!  spacing is computed.
!
!$OMP PARALLEL
      IF (.not.allocated(FinerDonor)) THEN
        allocate ( FinerDonor(Ngrids) )
        FinerDonor = 0
      END IF
!$OMP END PARALLEL
!
!  Logical switch indicating which coarser grid is a donor to a finer
!  grid external contact points. Their values are computed in
!  "metrics.F" after "CoarserDonor" is computed.
!
!$OMP PARALLEL
      IF (.not.allocated(DonorToFiner)) THEN
        allocate ( DonorToFiner(Ngrids) )
        DonorToFiner = .FALSE.
      END IF
!$OMP END PARALLEL
!
!  Number of refined grid time-steps. Their values are initialized to
!  the optimal values based on the spatial refinement ratio.  These
!  values are over-written in "metrics.F" to use the time-step
!  size ration between donor and receiver grid.  The user is allowed
!  to take larger divisible time-step with respect to the donor grid.
!  The user is responsible to set the appropriate refined grid time-step
!  for stability.
!
!$OMP PARALLEL
      IF (.not.allocated(RefineSteps)) THEN
        allocate ( RefineSteps(Ngrids) )
        DO ig=1,Ngrids
          RefineSteps(ig)=Refine_factor(ig)
        END DO
      END IF
!$OMP END PARALLEL
!
!  Refine grid time-steps counter with respect the coarse grid (ng=1)
!  single time-step.
!
      IF (.not.allocated(RefineStepsCounter)) THEN
        allocate ( RefineStepsCounter(Ngrids) )
        RefineStepsCounter = 0
      END IF
!
!  Counter for two-way exchange between fine and coarse grids.
!
      IF (.not.allocated(TwoWayCounter)) THEN
        allocate ( TwoWayCounter(Ngrids) )
        TwoWayCounter = 0
      END IF
!
!  Switch indicating which refined grid(s) include finer refined grids:
!  telescoping refinement.
!
!$OMP PARALLEL
      IF (.not.allocated(Telescoping)) THEN
        allocate ( Telescoping(Ngrids) )
        Telescoping = .FALSE.
      END IF
!$OMP END PARALLEL
!
!  Rolling index and time (seconds) used in the temporal interpolation
!  of contact point data.
!
!$OMP PARALLEL
      IF (.not.allocated(RollingIndex)) THEN
        allocate ( RollingIndex(Ncontact) )
        RollingIndex = 0
      END IF

      IF (.not.allocated(RollingTime)) THEN
        allocate ( RollingTime(2,Ncontact) )
        RollingTime = 0
      END IF
!$OMP END PARALLEL
!
!-----------------------------------------------------------------------
!  Allocate composite grids contact regions structure and arrays.
!-----------------------------------------------------------------------
!
      IF (ANY(CompositeGrid)) THEN
!
!  Allocate composite grid contact regions (type T_COMPOSITE) structure.
!
        allocate ( COMPOSITE(Ncontact) )
!
!  Allocate arrays in composite grids contact regions structure.
!
        DO cr=1,Ncontact
          dg=donor_grid(cr)

          allocate ( COMPOSITE(cr) % bustr(4,NpointsU(cr)) )
          allocate ( COMPOSITE(cr) % bvstr(4,NpointsV(cr)) )

          allocate ( COMPOSITE(cr) % ubar(4,NpointsU(cr),2) )
          allocate ( COMPOSITE(cr) % vbar(4,NpointsV(cr),2) )
          allocate ( COMPOSITE(cr) % zeta(4,NpointsR(cr),2) )

          allocate ( COMPOSITE(cr) % rzeta(4,NpointsR(cr)) )

# ifdef SOLVE3D
          allocate ( COMPOSITE(cr) % DU_avg1(4,NpointsU(cr)) )
          allocate ( COMPOSITE(cr) % DV_avg1(4,NpointsV(cr)) )
          allocate ( COMPOSITE(cr) % Zt_avg1(4,NpointsR(cr)) )

          allocate ( COMPOSITE(cr) % u(4,N(dg),NpointsU(cr)) )
          allocate ( COMPOSITE(cr) % v(4,N(dg),NpointsV(cr)) )

          allocate ( COMPOSITE(cr) % Huon(4,N(dg),NpointsU(cr)) )
          allocate ( COMPOSITE(cr) % Hvom(4,N(dg),NpointsV(cr)) )

          allocate ( COMPOSITE(cr) % t(4,N(dg),NpointsR(cr),NT(dg)) )
# endif
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Allocate refinement grids contact regions structure.
!-----------------------------------------------------------------------
!
      IF (ANY(RefinedGrid)) THEN
!
!  Allocate refinement grids contact region (type T_REFINED) structure.
!
        allocate ( REFINED(Ncontact) )
!
!  Allocate arrays in refinement grids contact region structure.
!
        DO cr=1,Ncontact
          rg=receiver_grid(cr)

          allocate ( REFINED(cr) % ubar(4,NpointsU(cr),2) )
          allocate ( REFINED(cr) % vbar(4,NpointsV(cr),2) )
          allocate ( REFINED(cr) % zeta(4,NpointsR(cr),2) )

          allocate ( REFINED(cr) % DU_avg2(4,NpointsU(cr),2) )
          allocate ( REFINED(cr) % DV_avg2(4,NpointsV(cr),2) )

          allocate ( REFINED(cr) % on_u(NpointsU(cr)) )
          allocate ( REFINED(cr) % om_v(NpointsV(cr)) )

# ifdef SOLVE3D
          allocate ( REFINED(cr) % u(4,N(rg),NpointsU(cr),2) )
          allocate ( REFINED(cr) % v(4,N(rg),NpointsV(cr),2) )

          allocate ( REFINED(cr) % t(4,N(rg),NpointsR(cr),2,NT(rg)) )
# endif
        END DO
      END IF

      RETURN
#else
      SUBROUTINE set_contact
#endif
      END SUBROUTINE set_contact
