.. _WRF_data:

================================
Running real cases with WRF data
================================

The CFBM can run in offline mode using WRF atmospheric fields as input. To run CFBM in the offline mode, users need to provide WRF atmospheric data, including all available timestamps, in a file named ```wrf.nc```, along with static inpupts in ```geo_em.d01.nc```. The atmospheric data must include wind components (U, V), geopotential heights (PH, PHB), surface variables used with the fuel moisture model (RAINC, RAINNC, T2, Q2, PSFC), and roughness length (ZNT). 

The simulation configuration is set up in the ```namelist.fire``` file, as described in Section 2.2. Users must specify the number of vertical levels (```kde```) and the time interval for incoming atmospheric data (```interval_atm```) to match the data provided in the ```wrf.nc``` file.

Users can obtain the model from the github repository:

```git clone https://github.com/NCAR/fire_behavior.git```

To compile the code on derecho, run:

```./compile.sh --env-auto```

If the compilation is successful, the model can be run using ```fire_behavior.exe``` located in the ```build``` directory.

An example is provided in ```test/test7/``` directory with the ```namelist.fire```.

If the simulation is successful, the model outputs are written to files named 'fire_output_*', and diagnostic messages are recorded in the 'log' file.

