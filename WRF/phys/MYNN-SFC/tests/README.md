## CI tests for MYNN-SFC

* Cases include:
	* ccpp_test(), on July 13, 2019 at SGP site over land
	* wrf_test('lnd'), on Aug 23, 2024 near (41.91, -73.12) over land
	* wrf_test('wat'), on Aug 23, 2024 near (41.78, -67.61) over water
	* wrf_test('icy'), on Feb 28, 2024 near (46.54, -73.53) over land covered with snow

* Input:

    * Input file for ccpp_test() is generated using CCPP SCM
    * Input files for wrf_test() is generated using WRFv4.5.1

* Output:
    * The CI test generates output files under ```tests/data``` dir

* To run the CI tests offline:

	```
	cd tests
	make
	./driver
	```