# MYNN-EDMF

The Mellor–Yamada–Nakanishi–Niino (MYNN) (Nakanishi and Niino 2001, 2004, 2006, and
2009) scheme has been adpted and developed (Olson et al. 2019) for use in NOAA's 
operational forecast models (RAP, HRRR, and RRFS). It has been integrated into several 
modeling framworks, such as the Advanced Research version of the Weather Research and 
Forecasting Model (WRF-ARW) (Skamarock et al. 2019), Common Community Physics Package 
(CCPP) (Bernardet et al. 2024), and Model Prediction Across Scales (MPAS) (Skamarock 
et al. 2012). To centralize the development of the MYNN-EDMF for all applications, this 
scheme code has been universalized and made accessible in this stand-alone submodule 
repository, which is then connected to each of the modeling frameworks mentioned above. 
All future development of the MYNN-EDMF that is intended for NOAA’s operational use 
will be hosted in this public-facing submodule repository.

The most up-to-date documentation of the MYNN-EDMF is new 2026 tech note:

Olson, Joseph B., Wayne M. Angevine, David D. Turner, Xia Sun, Julia M. Simonson,
Clark Evans, Jaymes S. Kenyon, Haiqin Li, Jordan Schnell, Franciano S. Puhales,
Tiziana Cherubini, Weiwei Li, and Man Zhang, 2026: A Description of the MYNN-EDMF
Turbulence Scheme. NOAA Tech. Memo. OAR GSL-77. 60 pp. doi:10.25923/rahr-sj70

The older tech note (Olson et al. 2019; doi:10.25923/n9wm-be49) still has some value,
for example when referencing older code in HRRRv4/RAPv5, but the updated tech note must
be used for referencing the contemporary code.

If you'd like to contribute, please reach out to the lead developer, Joseph Olson at
joseph.b.olson@noaa.gov, but keep in mind that all modification must be designed to work in
all 3 model frameworks (WRF, MPAS, and CCPP). Thanks in advance.
