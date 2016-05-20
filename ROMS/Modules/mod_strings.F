#include "cppdefs.h"
      MODULE mod_strings
!
!svn $Id: mod_strings.F 1368 2011-01-21 00:34:18Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  cdt         F90/F95 compiler used.                                  !
!  fflags      F90/F95 compiler flags.                                 !
!  title       Title of model run.                                     !
!  Coptions    Activated C-preprocessing options.                      !
!  Mregion     Start index of message passage code profiling region.   !
!  Nregion     Number of total code profiling regions.                 !
!                                                                      !
!  StateMsg    Model state processing messages:                        !
!                ( 1) Read state initial conditions                    !
!                ( 2) Read previous state initial conditions           !
!                ( 3) Read previous adjoint state solution             !
!                ( 4) Read latest adjoint state solution               !
!                ( 5) Read initial/model normalization factors         !
!                ( 6) Read correlation standard deviation              !
!                ( 7) Read impulse forcing                             !
!                ( 8) Read v-space increments                          !
!                ( 9) Read background state                            !
!                (10) Read boundary normalization factors              !
!                (11) Read forcing normalization factors               !
!                (12) Read surface forcing and/or OBS increments       !
!                                                                      !
!  Pregion     Model regions identifiers used for time profiling:      !
!                ( 1) Allocation and array initialization              !
!                ( 2) Ocean state initialization                       !
!                ( 3) Reading of input data                            !
!                ( 4) Processing of input data                         !
!                ( 5) Processing of output time averaged data          !
!                ( 6) Computation of vertical boundary conditions      !
!                ( 7) Computation of global information integrals      !
!                ( 8) Writing of output data                           !
!                ( 9) Model 2D kernel                                  !
!                (10) Lagrangian floats trajectories                   !
!                (11) Tidal forcing                                    !
!                (12) 2D/3D coupling, vertical metrics                 !
!                (13) Omega vertical velocity                          !
!                (14) Equation of state for seawater                   !
!                (15) Biological module, source/sink terms             !
!                (16) Sediment tranport module, source/sink terms      !
!                (17) Atmosphere-Ocean bulk flux parameterization      !
!                (18) KPP vertical mixing parameterization             !
!                (19) GLS vertical mixing parameterization             !
!                (20) My2.5 vertical mixing parameterization           !
!                (21) 3D equations right-side terms                    !
!                (22) 3D equations predictor step                      !
!                (23) Pressure gradient                                !
!                (24) Harmonic mixing of tracers, S-surfaces           !
!                (25) Harmonic mixing of tracers, geopotentials        !
!                (26) Harmonic mixing of tracers, isopycnals           !
!                (27) Biharmonic mixing of tracers, S-surfaces         !
!                (28) Biharmonic mixing of tracers, geopotentials      !
!                (29) Biharmonic mixing of tracers, isopycnals         !
!                (30) Harmonic stress tensor, S-surfaces               !
!                (31) Harmonic stress tensor, geopotentials            !
!                (32) Biharmonic stress tensor, S-surfaces             !
!                (33) Biharmonic stress tensor, geopotentials          !
!                (34) Corrector time-step for 3D momentum              !
!                (35) Corrector time-step for tracers                  !
!                (36) Two-way Atmosphere-Ocean models coupling         !
!                (37) Bottom boundary layer module                     !
!                (38) GST Analysis eigenproblem solution               !
!                (39) Multiple-grid nesting processing                 !
!                (40) Message Passage: 2D halo exchanges               !
!                (41) Message Passage: 3D halo exchanges               !
!                (42) Message Passage: 4D halo exchanges               !
!                (43) Message Passage: lateral boundary exchanges      !
!                (44) Message Passage: data broadcast                  !
!                (45) Message Passage: data reduction                  !
!                (46) Message Passage: data gathering                  !
!                (47) Message Passage: data scattering                 !
!                (48) Message Passage: boundary data gathering         !
!                (49) Message Passage: point data gathering            !
!                (50) Message Passage: multi-model coupling            !
!                                                                      !
!=======================================================================
!
        implicit none

        character (len=80)   :: title
        character (len=2048) :: Coptions

        integer, parameter :: Mregion = 40
        integer, parameter :: Iceregion = 51
        integer, parameter :: Nregion = 56

        character (len=55), dimension(12) :: StateMsg =                 &
     &    (/'Read state initial conditions,               ',            &
     &      'Read previous state initial conditions,      ',            &
     &      'Read previous adjoint state solution,        ',            &
     &      'Read latest adjoint state solution,          ',            &
     &      'Read initial/model normalization factors,    ',            &
     &      'Read correlation standard deviation,         ',            &
     &      'Read impulse forcing,                        ',            &
     &      'Read v-space increments,                     ',            &
     &      'Read background state,                       ',            &
     &      'Read boundary normalization factors,         ',            &
     &      'Read forcing normalization factors,          ',            &
     &      'Read surface forcing and/or OBS increments,  '/)

        character (len=50), dimension(Nregion) :: Pregion =             &
     &    (/'Allocation and array initialization ..............',       &
     &      'Ocean state initialization .......................',       &
     &      'Reading of input data ............................',       &
     &      'Processing of input data .........................',       &
     &      'Processing of output time averaged data ..........',       &
     &      'Computation of vertical boundary conditions ......',       &
     &      'Computation of global information integrals ......',       &
     &      'Writing of output data ...........................',       &
     &      'Model 2D kernel ..................................',       &
     &      'Lagrangian floats trajectories ...................',       &
     &      'Tidal forcing ....................................',       &
     &      '2D/3D coupling, vertical metrics .................',       &
     &      'Omega vertical velocity ..........................',       &
     &      'Equation of state for seawater ...................',       &
     &      'Biological module, source/sink terms .............',       &
     &      'Sediment tranport module, source/sink terms ......',       &
     &      'Atmosphere-Ocean bulk flux parameterization ......',       &
     &      'KPP vertical mixing parameterization .............',       &
     &      'GLS vertical mixing parameterization .............',       &
     &      'My2.5 vertical mixing parameterization ...........',       &
     &      '3D equations right-side terms ....................',       &
     &      '3D equations predictor step ......................',       &
     &      'Pressure gradient ................................',       &
     &      'Harmonic mixing of tracers, S-surfaces ...........',       &
     &      'Harmonic mixing of tracers, geopotentials ........',       &
     &      'Harmonic mixing of tracers, isopycnals ...........',       &
     &      'Biharmonic mixing of tracers, S-surfaces .........',       &
     &      'Biharmonic mixing of tracers, geopotentials ......',       &
     &      'Biharmonic mixing of tracers, isopycnals .........',       &
     &      'Harmonic stress tensor, S-surfaces ...............',       &
     &      'Harmonic stress tensor, geopotentials ............',       &
     &      'Biharmonic stress tensor, S-surfaces .............',       &
     &      'Biharmonic stress tensor, geopotentials ..........',       &
     &      'Corrector time-step for 3D momentum ..............',       &
     &      'Corrector time-step for tracers ..................',       &
     &      'Two-way Atmosphere-Ocean models coupling .........',       &
     &      'Bottom boundary layer module .....................',       &
     &      'GST Analysis eigenproblem solution ...............',       &
     &      'Multiple-grid nesting processing .................',       &
     &      'Message Passage: 2D halo exchanges ...............',       &
     &      'Message Passage: 3D halo exchanges ...............',       &
     &      'Message Passage: 4D halo exchanges ...............',       &
     &      'Message Passage: lateral boundary exchanges ......',       &
     &      'Message Passage: data broadcast ..................',       &
     &      'Message Passage: data reduction ..................',       &
     &      'Message Passage: data gathering ..................',       &
     &      'Message Passage: data scattering..................',       &
     &      'Message Passage: boundary data gathering .........',       &
     &      'Message Passage: point data gathering ............',       &
     &      'Message Passage: multi-model coupling ............',       &
     &      'Ice thermodynamics................................',       &
     &      'Ice rheology coefficients.........................',       &
     &      'Generate coefficients for ice dynamics solver.....',       &
     &      'Generate RHS for ice dynamics solver..............',       &
     &      'Iterative solver of ice dynamics..................',       &
     &      'Advection of ice tracers..........................'/)

        character (len=80) :: my_os = MY_OS
        character (len=80) :: my_cpu = MY_CPU
        character (len=80) :: my_fort = MY_FORT
        character (len=80) :: my_fc = MY_FC
        character (len=160) :: my_fflags = MY_FFLAGS

      END MODULE mod_strings
