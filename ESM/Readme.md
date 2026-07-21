
<img width="500" alt="image" src="https://github.com/myroms/roms/assets/23062912/2fa815f4-df51-4671-b9ed-188995b87a7b">


For detailed documentation and instructions, please check the following links:

| $\textcolor{red}{\textsf{Information}}$ | $\textcolor{red}{\textsf{WikiROMS}}$   |
|-------------------------------------|--------------------------------------------|
| General                             | https://www.myroms.org/wiki/Model_Coupling |
| **ROMS** Native **NUOPC** _cap_     | https://www.myroms.org/wiki/NUOPC_Cap      |
| **ROMS** Standalone **NUOPC** _cap_ | https://www.myroms.org/wiki/NUOPC_Cap_UFS  |


## Earth System Model (ESM) Coupling:

This directory contains several files used for multi-model coupling
using the Earth System Modeling Framework (**ESMF**) with the National
Unified Operational Prediction Capability (**NUOPC**) layer.

The **NUOPC** layer is a simplified **API** on top of the **ESMF** library
(version **8.0** or higher) that provides conventions and templates to facilitate
the smooth coupling between Earth System Model (**ESM**) components. **ROMS**
offers two distinct **NUOPC**-based modules: one for its native coupling framework
(**`esmf_roms.h`**) and the other as a standalone interface for third-party
coupling systems (**`cmeps_roms.h`**).  That is, **ROMS** coupling infrastructure
allows both **DRIVER** and **COMPONENT** methods of operation.

<p align="center">
    <img src="https://www.myroms.org/trac/ROMS_Coupling.png" width=50% height=50%>
</p>

In the **`DRIVER`** method, it provides all the interfaces needed to couple
to other **ESM** components, including the executable driver, **NUOPC**-based generic
**ESM** component services, model gridded components or **NUOPC** _cap_ modules,
connectors between components for re-gridding source and destination
fields, input scripts, and coupling metadata management.

A **NUOPC** model _cap_ is a Fortran code **API** layer that sits on top of the
**ESM** component, making calls to the numerical kernel via the _`initialize`_,
_`run`_, and _`finalize`_ computational phases.

Alternatively, in the **`COMPONENT`** method, the **NUOPC**-based **ROMS** _cap_
module **`Master/cmeps_roms.h`** is provided, and it can be adapted and
incorporated into other **NUOPC**-based coupling systems, like the
[UFS_coastal Framework](https://github.com/oceanmodeling/ufs-coastal)
using the Community Mediator for Earth Prediction Systems (**CMEPS**) and
the Community Data Models for Earth Predictive Systems (**CDEPS**) coupling
interfaces.

Currently, the following files are available for coupling with the **ESMF/NUOPC**
library:

| $\textcolor{blue}{\textsf{Coupling Files}}$ | $\textcolor{blue}{\textsf{Description}}$ |
|-----------------------|-------------|
| **cmeps_roms.h**      | **ROMS** standalone interface for the **UFS** using **CDEPS/CMEPS** |
| **coupler.F**         | **ESMF/NUOPC** or **MCT** native coupler module |
| **esmf_coupler.h**    | **ESMF** regridding operators using connectors |
| **esmf_atm.F**        | Atmosphere gridded model (**ATM**) component module |
| **esmf_atm_coamps.h** | **COAMPS** gridded component **NUOPC** layer module |
| **esmf_atm_regcm.h**  | **RegCM** gridded component **NUOPC** layer module |
| **esmf_atm_wrf.h**    | **WRF** gridded component **NUOPC** layer module |
| **esmf_data.F**       | Coupling **DATA** component used in incongruent grids coupling |
| **emsf_driver.h**     | Configures, Creates, Initialize, Run, and Finalizes the **ESM** coupling system |
| **esmf_esm.F**        | Sets **ESM** gridded components services shared-objects and RunSequence |
| **esmf_ice.F**        | Seaice gridded model (**SEA-ICE**) component module |
| **esmf_ice_cice6.h**  | **CICE6** gridded component **NUOPC** layer module |
| **esmf_roms.h**       | **ROMS** native coupling framework **NUOPC** layer module |
| **esmf_roms.F**       | **ROMS** native or standalone Ocean component module |
| **esmf_wav.F**        | Waves gridded component (**WAVE**) component module |
| **esmf_wav_wam.h**    | **WAM** gridded component **NUOPC** layer module |
| **mod_esmf_esm.F**    | **ROMS** native coupling framework support module |

## Coupling Design:

The strategy is to couple to other **ESM** components with no or minimal
changes to the code distributed by developers or repositories.  The User is
responsible for subscribing to those repositories, downloading, and installing
the other **ESM** components. The coupling with such **ESM** components is expected
not to be affected by its version (previous, current, or future) since the
**NUOPC**-based layer is generic.

However, sometimes, we need to circumvent technical problems when coupling
to other **ESM** components and provide build scripts to facilitate
compiling, linking, and correcting interference to deprecated libraries.
Therefore, this directory may also contain scripts and modified **ESM** component
files that substitute the ones distributed from source repositories
to solve such technical issues.


## WRF Coupling:

Previously, we needed to patch several **WRF** configuration files to correct the
NetCDF4 dependencies and rename collisions with newer versions of the **ESMF** library.
**WRF** adopted the **ESMF** date/time management **API** several years ago and kept the same
function names, but they were never updated. Nowadays, this patching is not required
since we loaded the changes to **feature/fix_compile_dependencies** branch of the
[WRF GitHub](https://github.com/wrf-model/WRF) repository, which started from
**WRF** version **4.5.1** on Oct. 10,  2023.

Use the following command to download **WRF** from GitHub:

```
git clone -b feature/fix_compile_dependencies https://github.com/wrf-model/WRF
```

The compiling of **WRF** is archaic and convoluted. Thus, we provide a build script
(**build_wrf.csh** or **build_wrf.sh**) to facilitate the compiling and linking in
various computer environments, including Spack-Stack. The resulting **WRF** libraries
can be either static or shared. Also, an additional script (**wrf_move.csh** or **wrf_move.sh**)
is supplied to move (**`-move`** option) **WRF** objects, libraries, and executables to the
User's Project sub-directory **`./Build_wrf`** to facilitate running **WRF** as an atmospheric
component in an **ESMF/NUOPC** coupling system and keep all the configurations in the same place.
To run **WRF**, we need an increasing number of formatted and unformatted input data files that
need to be located annoyingly in the Project directory. We use the **wrf_links.csh** or
**wrf_links.sh** to create the links to those files from the source code root directory.

For example, use the following command to compile and link **WRF**:

```
build_wrf.sh -j 10 -move
```

If compiling with **gfortran**, you may need to activate the environmental variable
**GFORTRAN_CONVERT_UNIT** with **big_endian** to avoid end-of-file when reading binary files
like **RRTMG_LW_DATA**:

```
  export GFORTRAN_CONVERT_UNIT='big_endian'
or
  setenv GFORTRAN_CONVERT_UNIT big_endian
```

## Files Summary:

We provide various scripts and files that can be used as an example and templates to
configure the **ESMF/NUOPC**-based coupling system:

| $\textcolor{blue}{\textsf{Coupling Files}}$ | $\textcolor{blue}{\textsf{Description}}$        |
|--------------------------------|--------------------------------------------------------------|
| **build_cice.csh**             | Script to compile and link **CICE6** gridded component       |
| **build_ufs.csh, .sh**         | **UFS** CMake compiling and linking CSH and BASH scripts     |
| **build_wps.csh, .sh**         | **WPS** compiling and linking CSH and BASH scripts           |
| **build_wrf.csh, .sh**         | **WRF** compiling and linking CSH and BASH scripts           |
| **wrf_links.chs, .sh**         | Creates links to **WRF** input data files in Project directory  |
| **wrf_move.csh, .sh**          | Moves **WRF** objects, libraries, and executables to Project directory |
|                                |                                                              |
| **coamps_explicit.runconfig**  | **DATA-COAMPS-ROMS** ESMF RunSequence explicit configuration |
| **coamps_implicit.runconfig**  | **DATA-COAMPS-ROMS** ESMF RunSequence implicit configuration |
| **data.runconfig**             | **DATA-ROMS** ESMF RunSequence configuration                 |
| **data_snapshots.runconfig**   | **DATA-ROMS** ESMF RunSequence configuration with snapshots  |
|                                |                                                              |
| **coupling_esmf.in**           | Native framework input script: **`mpirun -np 8 romsM coupling_esmf.in > & log &`** |
| **coupling_esmf.yaml**         | Generic native framework exchange fields metadata            |
| **coupling_esmf_coamps.yaml**  | **DATA-COAMPS_ROMS** system exchange fields metadata         |
| **coupling_esmf_wrf.yaml**     | **DATA-WRF-ROMS** system exchange fields metadata            |
| **cmeps_roms.yaml**            | **UFS-ROMS** coupling configuration with **CDEPS/CMEPS**     |

## ROMS-UFS Applications:

The BASH script **roms_test_ufs.sh** is provided to **`checkout`** just the directories
needed for configuring and running **ROMS-UFS** test applications from the **roms_test**
repository. Currently, only the [Hurricane Irene](https://github.com/myroms/roms_test/tree/main/IRENE)
application is available. Also, it illustrates the use of the **sparse-checkout** feature of git.
