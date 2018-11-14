# README for the COAWST Projects/Inlet_test/Coupled case

## Mark Hadfield, 2018-02-13

This case requires only SWAN and ROMS so can be built with the roms_build
script:

https://github.com/hadfieldnz/roms-scripts-mgh/blob/master/roms_build

The simulation-specific environment variables required by roms_build are
specified in the file build_parameters.

To build and run the case, set up host-specific variables like FORT and the
NETCDF variables, and invoke coawst_run, which calls roms_build to build the
executable and then runs

The coawst_run script accepts 2 options:

- -g sets the ROMS USE_DEBUG variable, which activates things like bounds
  checking and traceback. With this option selected, execution will be very
  slow.
  
- -s causes a SLURM script to be created and submitted for execution.
