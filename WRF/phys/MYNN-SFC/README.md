# MYNN-SFC

The Mellor–Yamada–Nakanishi–Niino (MYNN) surface layer scheme has been developed
for use in NOAA's operational forecast models (RAP, HRRR, and RRFS). It has been
integrated into several modeling framworks, such as the Advanced Research version of
the Weather Research and Forecasting Model (WRF-ARW) (Skamarock et al. 2019), Common
Community Physics Package (CCPP) (Bernardet et al. 2024), and Model Prediction Across
Scales (MPAS) (Skamarock et al. 2012). To centralize the development of the MYNN
surface layer scheme for all applications, this scheme code is currently being universalized
and made accessible in this stand-alone submodule repository, which is then connected to
each of the modeling frameworks mentioned above. All future development of the MYNN
surface layer scheme that is intended for NOAA’s operational use will be hosted in this
public-facing submodule repository.

The most up-to-date documentation of the MYNN surface layer scheme is:
https://repository.library.noaa.gov/view/noaa/30605
but many new features are being consolidated into this single universalized submodule,
so an updated tech note is targeted for publication by end of 2026.

If you'd like to contribute, please reach out to the lead developers:
Joseph Olson: joseph.b.olson@noaa.gov
Xia Sun: xia.sun@noaa.gov
