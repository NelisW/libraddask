# librad - An interface to libRadtran
## Introduction
libRadtran (http://www.libradtran.org/) is a collection of tools for atmospheric radiative transfer computation. 
The main driver utility for libRadtran is an executable called uvspec. This utility
reads an input file (generally *.INP) of keywords and keyword parameters that set up
and execute the radiative transfer computation. A variety of radiative transfer solvers
are provided with libRadtran, including various implementations of the classical Discrete Ordinates
Radiative Transfer (DISORT) method.

# MORTICIA and libRadtran
For MORTICIA, libRadtran is used to calculate radiances and irradiances for establishing the
radiative environment of a target and the amount of light emitted (longwave/thermal bands) or reflected
(shortwave/solar bands) from the target that reaches the sensor.

For rapid shortwave spectrum work, the band parametrisation by Kato is typically used in MORTICIA.
The Kato model has only 32 spectral bands across the whole available spectrum and all computations
are reduced to a vector of values representative of the Kato bands.

For more accurate work, such as simulations of satellite views of water targets, the REPTRAN
model is preferred. REPTRAN provides a full spectral calculation or rapid band model calculations
for a variety of satellite sensors.

To execute full simulations using MORTICIA it is necessary to have a working installation of
libRadtran.

# Twilight
libRadtran is able to calculate radiances and irradiances when the sun is below the horizon by up
to 9 degrees or more. This is well into "nautical" twilight. The MYSTIC monte carlo solver is
required to perform twilight computations.

# libRadtran Downloads
Further information and downloads are available at [libRadtran website](http://www.libradtran.org).

# radute - Radiometry Utilities
The moritica.rad package also inlcudes a number of utilities related to radiometry. This includes

- Creating, reading and writing of MODTRAN-style .flt files for specification of sensor Spectral Response
    functions (SRFs)
- Conversion of .flt class objects to lists of xarray.DataArray objects for use in the Spectra class.
- The RadEnv class for computation and handling of Radiant Environment Maps (REMs).





 

 
 
 

