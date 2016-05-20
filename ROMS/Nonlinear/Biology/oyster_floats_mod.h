!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Diego A. Narvaez   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for oyster model (Dekshenieks et al. 1993, 1996, 1997;   !
!  Narvaez et al. 2012a,b):                                            !
!                                                                      !
!  Larvae_GR0      growth rate (um/day)                                !
!  Larvae_size0    size in terms of length (um)                        !
!  food_supply     food supply (constant source; mg Carbon/l)          !
!  settle_size     settlement size (um)                                !
!                                                                      !
!  Turbidity effect parameters on planktonic larvae growth:            !
!                                                                      !
!  turb_ambi       ambient turbidity level (g/l)                       !
!  turb_axis       turbidity linear axis crossing (g/l)                !
!  turb_base       turbidity exponential base factor (g/l)             !
!  turb_crit       critical turbidity value (g/l)                      !
!  turb_mean       turbidity exponential mean (g/l)                    !
!  turb_rate       turbidity exponential rate (1/(g/l))                !
!  turb_size       minimum larvae size (um) affected by turbidity      !
!  turb_slop       turbidity linear slope (1/(g/l))                    !
!                                                                      !
!  Planktonic larvae vertical migration (swimming) parameters:         !
!                                                                      !
!  slope_Sdec      coefficient due to decreasing salinity              !
!  slope_Sinc      coefficient due to increasing salinity              !
!  swim_Sdec       fraction active due to decreasing salinity          !
!  swim_Sinc       fraction active due to increasing salinity          !
!  swim_Tmax       maximum swimming time fraction                      !
!  swim_Tmin       minimum swimming time fraction                      !
!                                                                      !
!  Planktonic larvae sinking parameters:                               !
!                                                                      !
!  sink_base       exponential base factor (mm/s)                      !
!  sink_rate       exponential rate factor (1/um)                      !
!  sink_size       exponential mean size (um)                          !
!                                                                      !
!  Planktonic larvae swimming speed (mm/s) as a function of larval     !
!  size (um) and temperature (Celsius), look table parameters:         !
!                                                                      !
!  swim_DL         larval size J-axis increment                        !
!  swim_DT         temperature I-axis increment                        !
!  swim_Im         number of values in temperature I-axis              !
!  swim_Jm         number of values in larval size J-axis              !
!  swim_L0         starting value for larval size J-axis               !
!  swim_T0         starting value for temperature I-axis               !
!  swim_table      larval swimming speed look table                    !
!                                                                      !
!  Planktonic larvae growth rate factor (nondimensional) as a function !
!  salinity and temperature, look table parameters:                    !
!                                                                      !
!  Gfactor_DS      salinity I-axis increment                           !
!  Gfactor_DT      temperature J-axis increment                        !
!  Gfactor_Im      number of values in salinity I-axis                 !
!  Gfactor_Jm      number of values in temperature J-axis              !
!  Gfactor_S0      starting value for salinity I-axis                  !
!  Gfactor_T0      starting value for temperature J-axis               !
!  Gfactor_table   larval growth rate factor look table                !
!                                                                      !
!  Planktonic larvae growth rate as a function of food supply and      !
!  larval size, look table parameters:                                 !
!                                                                      !
!  Grate_DF        food supply I-axis increment                        !
!  Grate_DL        larval size J-axis increment                        !
!  Grate_Im        number of values in food supply I-axis              !
!  Grate_Jm        number of values in larval size J-axis              !
!  Grate_F0        starting value for food supply I-axis               !
!  Grate_L0        starting value for larval size J-axis               !
!  Grate_table     larval growth rate look table                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Oyster floats model parameters.
!
      integer  :: Gfactor_Im
      integer  :: Gfactor_Jm

      integer  :: Grate_Im
      integer  :: Grate_Jm

      integer  :: swim_Im
      integer  :: swim_Jm
!
      real(r8) :: Gfactor_DS
      real(r8) :: Gfactor_DT
      real(r8) :: Gfactor_S0
      real(r8) :: Gfactor_T0

      real(r8) :: Grate_F0
      real(r8) :: Grate_L0
      real(r8) :: Grate_DF
      real(r8) :: Grate_DL

      real(r8) :: swim_DL
      real(r8) :: swim_DT
      real(r8) :: swim_L0
      real(r8) :: swim_T0
!
      real(r8), allocatable :: Larvae_GR0(:)        ! um/day
      real(r8), allocatable :: Larvae_size0(:)      ! um
      real(r8), allocatable :: food_supply(:)       ! mg Carbon/l
      real(r8), allocatable :: settle_size(:)       ! um
!
      real(r8), allocatable :: sink_base(:)         ! mm/s
      real(r8), allocatable :: sink_rate(:)         ! 1/um
      real(r8), allocatable :: sink_size(:)         ! um
!
      real(r8), allocatable :: slope_Sinc(:)
      real(r8), allocatable :: slope_Sdec(:)
      real(r8), allocatable :: swim_Sdec(:)
      real(r8), allocatable :: swim_Sinc(:)
      real(r8), allocatable :: swim_Tmax(:)
      real(r8), allocatable :: swim_Tmin(:)
!
      real(r8), allocatable :: turb_ambi(:)         ! g/l
      real(r8), allocatable :: turb_axis(:)         ! g/l
      real(r8), allocatable :: turb_base(:)         ! g/l
      real(r8), allocatable :: turb_crit(:)         ! g/l
      real(r8), allocatable :: turb_mean(:)         ! g/l
      real(r8), allocatable :: turb_rate(:)         ! l/g
      real(r8), allocatable :: turb_size(:)         ! um
      real(r8), allocatable :: turb_slop(:)         ! l/g
!
      real(r8), allocatable :: Gfactor_table(:,:)
      real(r8), allocatable :: Grate_table(:,:)
      real(r8), allocatable :: swim_table(:,:)
!
      CONTAINS
!
      SUBROUTINE allocate_behavior
!
!=======================================================================
!                                                                      !
!  This routine allocates oyster model variables.                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(Larvae_GR0)) THEN
        allocate ( Larvae_GR0(Ngrids) )
      END IF

      IF (.not.allocated(Larvae_size0)) THEN
        allocate ( Larvae_size0(Ngrids) )
      END IF

      IF (.not.allocated(food_supply)) THEN
        allocate ( food_supply(Ngrids) )
      END IF

      IF (.not.allocated(settle_size)) THEN
        allocate ( settle_size(Ngrids) )
      END IF

      IF (.not.allocated(sink_base)) THEN
        allocate ( sink_base(Ngrids) )
      END IF

      IF (.not.allocated(sink_rate)) THEN
        allocate ( sink_rate(Ngrids) )
      END IF

      IF (.not.allocated(sink_size)) THEN
        allocate ( sink_size(Ngrids) )
      END IF

      IF (.not.allocated(slope_Sdec)) THEN
        allocate ( slope_Sdec(Ngrids) )
      END IF

      IF (.not.allocated(slope_Sinc)) THEN
        allocate ( slope_Sinc(Ngrids) )
      END IF

      IF (.not.allocated(swim_Sdec)) THEN
        allocate ( swim_Sdec(Ngrids) )
      END IF

      IF (.not.allocated(swim_Sinc)) THEN
        allocate ( swim_Sinc(Ngrids) )
      END IF

      IF (.not.allocated(swim_Tmax)) THEN
        allocate ( swim_Tmax(Ngrids) )
      END IF

      IF (.not.allocated(swim_Tmin)) THEN
        allocate ( swim_Tmin(Ngrids) )
      END IF

      IF (.not.allocated(turb_ambi)) THEN
        allocate ( turb_ambi(Ngrids) )
      END IF

      IF (.not.allocated(turb_axis)) THEN
        allocate ( turb_axis(Ngrids) )
      END IF

      IF (.not.allocated(turb_base)) THEN
        allocate ( turb_base(Ngrids) )
      END IF

      IF (.not.allocated(turb_crit)) THEN
        allocate ( turb_crit(Ngrids) )
      END IF

      IF (.not.allocated(turb_mean)) THEN
        allocate ( turb_mean(Ngrids) )
      END IF

      IF (.not.allocated(turb_rate)) THEN
        allocate ( turb_rate(Ngrids) )
      END IF

      IF (.not.allocated(turb_size)) THEN
        allocate ( turb_size(Ngrids) )
      END IF

      IF (.not.allocated(turb_slop)) THEN
        allocate ( turb_slop(Ngrids) )
      END IF

      RETURN
      END SUBROUTINE allocate_behavior
