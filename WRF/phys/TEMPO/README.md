TEMPO-v3.0.4

## Thompson-Eidhammer Microphysics Parameterization for Operations

[TEMPO documentation](https://ncar.github.io/TEMPO/)

### About
This repository houses the Thompson-Eidhammer Microphysics Parameterization for Operations, or TEMPO. TEMPO predicts clouds and precipitation in Earth-system models and is a continuation of the Thompson microphysics parameterization that is a popular choice in the Weather Research and Forecasting (WRF) model and run operationally in the High-Resolution Rapid Refresh (HRRR). TEMPO has been designed as a submodule so that development and changes are saved in a commit history, independent of a specific modeling framework.

The `main` branch of TEMPO is running in RRFSv2 development systems that use a UFS-community fork of the Model for Prediction Across Scales (MPAS[^1]) as a dynamical core.

Several key publications describe the research and model development that is contained in TEMPO. Note that code changes do occur -- mathematical formulations may be different in the code than as described in publications. When using TEMPO, some or all of these publications should be referenced. Additionally, commit hashes or tags should be included in publications.

**Thompson microphysics**
- Thompson, G., R. M. Rasmussen, and K. Manning, 2004: Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part I: Description and sensitivity analysis. Mon. Wea. Rev., 132, 519–542, https://doi.org/10.1175/1520-0493(2004)132<0519:EFOWPU>2.0.CO;2.
- Thompson, G., P. R. Field, R. M. Rasmussen, and W. D. Hall, 2008: Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part II: Implementation of a new snow parameterization. Mon. Wea. Rev., 136, 5095–5115, https://doi.org/10.1175/2008MWR2387.1.
- Thompson, G., 2013: High-resolution winter simulations of winter precipitation over the Colorado Rockies. Workshop on Parametrization of Clouds and Precipitation, Shinfield Park, Reading, ECMWF, 35–46, https://www.ecmwf.int/node/12672.
  
**Aerosol-aware microphysics**
- Thompson, G., and T. Eidhammer, 2014: A study of aerosol impacts on clouds and precipitation development in a large winter cyclone. J. Atmos. Sci., 71, 3636–3658, https://doi.org/10.1175/JAS-D-13-0305.1.

**2-moment graupel with predicted density (hail-aware) microphysics**
- Jensen, A. A., G. Thompson, K. Ikeda, and S. A. Tessendorf, 2023: Improving the Representation of Hail in the Thompson Microphysics Scheme. Mon. Wea. Rev., 151, 2307–2332, https://doi.org/10.1175/MWR-D-21-0319.1.

[^1]: https://github.com/ufs-community/MPAS-Model

### Working Practices and Governance
See the [Wiki](https://github.com/NCAR/TEMPO/wiki/Governance) for working practices and code governance.
