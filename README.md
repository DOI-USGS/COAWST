README
==========

This software is COAWST version 3.9. The user is recommended to read the User Manual for a full description of the procedures for installation, compilation, running, and visualization of the model. There have been several Trainings for the model, and instructions for downloading the tutorials are included in the manual.

Major updates for v3.9 include:
- Component updates to MCT v2.11, ROMS v4.3, WRF v4.8.0, WPS v4.7, SWAN v41.51A, WW3 v 7.14.2, InWave v 2.0;
- Wave number dependent wav-current interaction based on Banihashemi 2017.
- Coupling of multiple nests on the same refined level (WRF, ROMS, WW3, and SWAN).
- SCRIP automatic generation of weights file.
- Wave dissipation computed in component x- and y- directions.
- Z0_WAV computes roughness from wav model to atm based on the integral of the wave Source input (Sin) term.
- Z0_PORCHETTA computes roughness from wave model to atm based on wave directions.
- WAV2OCN fluxes computes surface stress sent to ocn model based on wave Source input (Sin) term. Provides complete separation of wave stresses into surface stress, depth-limited breaking, whitecapping, and bottom friction.
- InWave transform metrics for curvgrids.
- FLOATS WINDS and FLOATS_STOKES add wind slippage and stokes velocities to Eulerian drifter currents.
- Updated information in User manual. Please use updated inputs files, build scripts, and Compilers files.


# [COAWST Wiki](https://github.com/DOI-USGS/COAWST/wiki)
- [COAWST Trainings](https://github.com/DOI-USGS/COAWST/wiki/COAWST-Trainings)
- [Publications](https://github.com/DOI-USGS/COAWST/wiki/Publications)

