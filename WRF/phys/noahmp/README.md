# noahmp
Noah-MP Community Repository


This is the official Noah-MP unified Github repository for code downloading and contribution. Note that Noah-MP is a community model contributed by the whole Noah-MP scientific community. For maintenance and release of this GitHub, please contact: Cenlin He (cenlinhe@ucar.edu) and Fei Chen (feichen@ucar.edu).

Some changes have been made to the structure of archiving the stand-alone version of Noah-MP/HRLDAS code in the Github repository. Now, it separately archives the core Noah-MP source code, Noah-MP driver, and parameter tables in this unified Noah-MP Github repository and the rest of the HRLDAS related files (e.g., module_sf_urban.F, etc.) in another unified HRLDAS Github repository (https://github.com/NCAR/hrldas). The HRLDAS Github repo is already linked to this unified core Noah-MP code repository. This new archiving structure will allow different host model platforms/systems (e.g., HRLDAS, WRF, UFS, NWM, LIS, etc.) to connect to the core Noah-MP source code and develop their own host model drivers. 

Model developers can make code development based on the develop branch and create pull request to the develop branch. The pull request will be reviewed by Noah-MP model release team and if everything looks good, the new code development will be merged to the develop branch. Eventually, the updates in the develop branch will be merged to the master branch for official Noah-MP model release.

Branch structures of this Noah-MP repository:

1. "master" branch: most stable & recent version, updated whenever there are bug fixes or major model update/release;

2. "develop" branch: used for continuous NoahMP development, keep being updated by including bug fixes and code updates (e.g., new physics options, processes, etc.); 

3. "develop_no_gecros" branch: same as the "develop" branch, except for excluding the gecros crop section (to be consistent with recent Noah-MP changes for National Water Model (NWM));

4. "release-v4.3-WRF" branch: version archive of the stable version consistent with the WRFv4.3 release;

5. "release-v4.3-NWM" branch: version archive, same as the "release-v4.3-WRF" branch, except for excluding the gecros crop section (to be consistent with recent Noah-MP changes for National Water Model (NWM));

Some suggestions for model developers to contribute to Noah-MP code through the Github repository (typical procedures):

1. Step (1) Create a fork of this official NoahMP repository to your own Github account; 

2. Step (2) Make code updates based on the "develop" branch of the forked repository under your own account; 

3. Step (3) Finalize and test the code updates you make; 

4. Step (4) Submit a pull request for your model updates from your own forked Github repository to the "develop" branch of the official repository;

5. Step (5) The Noah-MP release team reviews and tests the model updates in the submitted pull request and discusses with the developer if there is any problem; 

6. Step (6) The Noah-MP release team confirms the pull request and merges the updated code to the "develop" and "develop_no_gecros" branches in the official repository;

7. Step (7) The Noah-MP release team will merge the updated "develop" branch to the master branch and version-release branch during the annual model release.

Note: If updates are made to both the core NoahMP source codes and other HRLDAS files, two separate pull requests need to be submitted to this NoahMP repository and the HRLDAS repository (https://github.com/NCAR/hrldas), respectively, regarding the changes in the Noah-MP code files and other HRLDAS-related files. This could be done by using the same titles for the two pull requests to ensure that the submitted code changes are handled together by the release team in the two repositories.
