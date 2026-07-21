# ROMS Metadata

# License

**Copyright (c) 2002-2026 The ROMS Group**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Overview

This directory contains several **CDL** files showing **ROMS** input **NetCDF**
file structure. There is a lot of software out there to generate
such input files. It is challenging to write generic programs
because it depends on the application. However, there is a straightforward
way to generate these files using the **CDL** and the **NetCDF** **`ncgen`**
program.

The **ROMS** metadata design is very rich and extensive. See "varinfo.yaml"
for a list of all the variables' names, units, attributes, associated time
variables, and scale factors. This is a user-friendly file, and variable
parameters can be quickly changed. Some users like
to change the long-name attribute to a language other than English to
facilitate automatic labeling during plotting. However, for portability
I think you should use the provided field variable name.

Currently, you can find the following **CDL** scripts:

``` 
    grd_spherical.cdl        Spherical grid NetCDF file

    ini_hydro.cdl            Initial conditions NetCDF file (hydrodynamics)
    ini_fennel.cdl           Initial conditions NetCDF file (hydrodynamics and biology)
    ini_ecosim.cdl           Initial conditions NetCDF file (hydrodynamics and bio-optics)
    ini_sed.cdl              Initial conditions NetCDF file (hydrodynamics and sediment)

    clm_ts.cdl               Temperature-Salinity climatology NetCDF file

    frc_uvstress.cdl         Forcing NetCDF file (surface momentum stresses)
    frc_fluxclm.cdl          Forcing NetCDF file (climatological heat fluxes variables)
    frc_bulk.cdl             Forcing NetCDF file (atmospheric variable for bulk fluxes)

    frc_rivers.cdl           Forcing NetCDF file (River point/sources)
    frc_tides.cdl            Forcing NetCDF file (tidal elevation and currents)

    bry_limit.cdl            Boundary NetCDF file (various time dimensions)
    bry_unlimit.cdl          Boundary NetCDF file (unlimited time dimensions)

    adsend.cdf               Adjoint sensitivity functional

    s4dvar_obs.cdl           4D-Var observations

    s4dvar_std_m.cdl         4D-Var model error covariance standard deviation
    s4dvar_std_i.cdl         4D-Var initial conditions error covariance standard deviation
    s4dvar_std_b.cdl         4D-Var open boundaries error covariance standard deviation
    s4dvar_std_f.cdl         4D-Var surface forcing error covariance standard deviation
```

Currently, there are two vertical, terrain-following coordinates
transformations in **ROMS**.  You need to choose the appropriate
**`standard_name`** attribute:

- Original transformation: **`ocean_s_coordinate_g1`**

``` nc
        double s_rho(s_rho) ;
                s_rho:long_name = "S-coordinate at RHO-points" ;
                s_rho:valid_min = -1. ;
                s_rho:valid_max = 0. ;
                s_rho:positive = "up" ;
                s_rho:standard_name = "ocean_s_coordinate_g1" ;
                s_rho:formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" ;

        double s_w(s_w) ;
                s_w:long_name = "S-coordinate at W-points" ;
                s_w:valid_min = -1. ;
                s_w:valid_max = 0. ;
                s_w:positive = "up" ;
                s_w:standard_name = "ocean_s_coordinate_g1" ;
                s_w:formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc" ;
```

- New transformation: **`ocean_s_coordinate_g2`**

``` nc
        double s_rho(s_rho) ;
                s_rho:long_name = "S-coordinate at RHO-points" ;
                s_rho:valid_min = -1. ;
                s_rho:valid_max = 0. ;
                s_rho:positive = "up" ;
                s_rho:standard_name = "ocean_s_coordinate_g2" ;
                s_rho:formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" ;

        double s_w(s_w) ;
                s_w:long_name = "S-coordinate at W-points" ;
                s_w:valid_min = -1. ;
                s_w:valid_max = 0. ;
                s_w:positive = "up" ;
                s_w:standard_name = "ocean_s_coordinate_g2" ;
                s_w:formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc" ;
```
Notice that the nondimensional, fractional, stretched vertical coordinate
except for the attribute's value **`standard_name`** is the same.

You can easily edit any of these files to change the **NetCDF** file name, change
dimensions, add and remove variables, and add and modify global attributes.
A **NetCDF** file can be created by typing:

``` d
    ncgen -b my_file.cdl
```

Then, you can use any program to write your data into the created **NetCDF**
file. 

Notice that **ROMS** allows multiple forcing **NetCDF** files. See input script.


