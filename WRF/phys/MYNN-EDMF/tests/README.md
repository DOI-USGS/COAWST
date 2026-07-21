## CI tests for MYNN-EDMF
* Cases include:
    * clr, land point (43.40°N, 73.02°W) on Aug 23, 2024 
    * non-precip StCu, ocean point (40.84°N, 66.80°W) on Mar 07, 2024
* Input:

    * Input files are generated using WRF

* To run the CI tests offline:
    ```
    export NFHOME=/path/to/your/netcdf/dir
    cd tests
    make
    ./driver
    ```
