.. _SRW:

================================================
Coupling to the UFS: Running the CFBM in the SRW
================================================

This schematic figure illustrates the integration of the Community Fire Behavior Model (CFBM) as a component within the Unified Forecast System (UFS) using the NUOPC (National Unified Operational Prediction Capability) framework.

The blue box at the top represents the UFS atmosphere component, which provides atmospheric variables (e.g., wind, temperature) to drive the fire model. The orange box represents CFBM that simulates fire spread and calculates heat, moisture, and fire smoke, which feed back into the atmospheric component. The NUOPC cap acts as an interface, facilitating the exchange of data between the UFS atmosphere and the CFBM. The CFBM requires static input in the form of geo_em.d01, which contains refined grids, fire fuels, and detailed topography to accurately simulate fire spread.

.. !.. image:: https://github.com/NCAR/fire_behavior/blob/develop/doc/CFBM-NUOPC.jpeg
.. !  :width: 400
.. !  :alt: CFBM-NUOPC coupling schematic
.. !  :align: center

UFS Short-Range Weather Application (SRW)
=========================================

The CFBM has been coupled to the UFS Weather Model for both one-way (atmosphere -> fire) and two-way coupled simulations.
Simulations can be run using the UFS Short-Range Weather Application (SRW), a community-supported application for running
numerical weather prediction simulations on limited-area domains. For information on using this capability, see
`the SRW Users Guide <https://ufs-srweather-app.readthedocs.io/en/latest/UsersGuide/BuildingRunningTesting/FIRE.html>`_.

A preprint on scientific results using this new UFS Fire capability is available (:cite:`CFBM`).

