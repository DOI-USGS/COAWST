#include "cppdefs.h"
#if defined FOUR_DVAR
!!
!!git $Id$
!!================================================= Hernan G. Arango ===
!!  Copyright (c) 2002-2026 The ROMS Group                             !
!!    Licensed under a MIT/X style license                             !
!!    See License_ROMS.md                                              !
!!======================================================================
!!                                                                     !
!! This module computes normalization factors to model the spreading   !
!! of the background-error symmetric covariance matrix (B). These      !
!! coefficients ensure that the diagonal elements of B are equal to    !
!! unity. In applications involving land or sea masking, the           !
!! normalization factors can significantly modify the covariance       !
!! structure near boundaries.                                          !
!!                                                                     !
!! Spreading is modeled using pseudo-diffusion operators in spatial    !
!! correlation space, where diffusion coefficients (K) are proportional!
!! to the square of the correlation length scale (Daley, 1992).        !
!! Correlation scales may be constant or spatially varying for each    !
!! variable in the control vector. The background-error covariance     !
!! matrix (B) can be represented using either a multiscale or monoscale!
!! approach. In the multiscale formulation, B is expressed as a linear !
!! combination of distinct spatial scales, ranging from large to small.!
!! This approach enables the representation of both broad structures   !
!! and fine-scale features, while reducing scale aliasing in the data  !
!! assimilation cost function (Weaver et al., 2016). In contrast, the  !
!! monoscale formulation employs a single scale.                       !
!!                                                                     !
!! The multiscale approach uses an implicit horizontal pseudo-         !
!! diffusion operator implemented via Conjugate Gradient (CG) and      !
!! Chebyshev Iterations (CI). Conversely, the monoscale default method !
!! employs an explicit horizontal algorithm. Error correlations are    !
!! considered separable in the horizontal and vertical directions.     !
!! The vertical diffusion operator remains implicit in both the        !
!! monoscale and multiscale approaches.                                !
!!                                                                     !
!! The multiscale concept, expressed as B = SUM(Wi*Bi), represents a   !
!! weighted sum of B values corresponding to various scales. Typically,!
!! two to four different scales are combined in practical applications.!
!! The weight coefficients Wi are constrained to sum to unity. The     !
!! resulting correlation functions belong to the Matern-class family,  !
!! which accommodates complex functional shapes.                       !
!!                                                                     !
!! Horizontally spatially varying correlation length scales for the    !
!! K-diffusion coefficient are permitted. A horizontal map of isotropic!
!! or anisotropic correlations can be specified and is read from an    !
!! input NetCDF file. For computational efficiency, these correlations !
!! require an implicit diffusion CG/CI solver.                         !
!!                                                                     !
!! References:                                                         !
!!                                                                     !
!!   Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!!     sphere using a generalized diffusion equation, Q.J.R. Meteorol. !
!!     Soc, 127, 1815-1846, doi:10.1002/qj.49712757518.                !
!!                                                                     !
!!   Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation  !
!!     operators based on an implicitly formulated diffusion equation  !
!!     solved with the Chebyshev iteration, Q.J.R. Meteorol. Soc.,     !
!!     142, 455-471, doi:10.1002/qj.2664.                              !
!!                                                                     !
!!======================================================================
!!
# ifdef MULTI_SCALE_B
#  include "normalization_multi.h"
# else
#  include "normalization_mono.h"
# endif
#else
      MODULE normalization_mod
      END MODULE normalization_mod
#endif
