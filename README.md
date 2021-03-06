# libraddask


## Introduction
libRadtran (http://www.libradtran.org/) is a collection of tools for atmospheric radiative transfer computation. 
The main driver utility for libRadtran is an executable called uvspec. This utility
reads an input file (generally *.INP) of keywords and keyword parameters that set up
and execute the radiative transfer computation. A variety of radiative transfer solvers
are provided with libRadtran, including various implementations of the classical Discrete Ordinates
Radiative Transfer (DISORT) method.
Further information and downloads are available at [libRadtran website](http://www.libradtran.org).

`libraddask` comprises three main parts:

1. `libraddask.rad` includes copies of the libRadtran files contained in `libradtran-2.0.3/src_py`, but converted to Python 3 status.

2. `libraddask.rad` also includes a number of utilities related to radiometry.  This code was originally written by Derek Griffith (https://github.com/derekjgriffith/MORTICIA) for Python 2. The Python 3 converted files are included here. These include

    - Creating, reading and writing of MODTRAN-style .flt files for specification of sensor Spectral Response functions (SRFs)
    - Conversion of .flt class objects to lists of xarray.DataArray objects for use in the Spectra class.
    - The RadEnv class for computation and handling of Radiant Environment Maps (REMs).


3. Documentation and examples to use Dask for distributed libRadTran runs (i.e., on remote clusters).


## Tools

The Python-3-converted `libRadtran/src_py` Python files are in the  `libraddask\libradtran` folder.  At the time of writing (20200420) the files in this folder are based on the libRadtran version 2.0.2 Python files, the port of the 2.0.3 files are still in process.

The utility files are in the `rad/` folder.


## Documentation
The documentation is provided in LaTeX format in the `doc/` folder.
The PDF version of the documentation might not be built in early beta releases, please build it yourself.

The `doc/` folder also contains a Jupyter notebook that demonstrates the use of the `libraddask` package.

## Installation and Requirements

`libraddask` requires Python 3.

A working installation of [libRadtran](http://www.libradtran.org) is required to compute radiant environment  maps and atmospheric transmittance. In the Monte carlo/statistical mode of operation, a compute cluster is generally required to achieve reasonable execution time of large numbers of runs.

Parallel computation is performed  using the [`Dask`](https://docs.dask.org/en/latest/) package, which works in  IPython/Jupyter notebooks as well as other Python launch modes.
 
`libraddask` has been developed and tested using the [Anaconda](https://www.continuum.io/downloads) distribution from [Continuum Analytics](https://www.continuum.io/).
 
## Repository
 The master repository for `libraddask` is  hosted on [GitHub](http://www.github.org) at 
 [https://github.com/NelisW/libraddask](https://github.com/NelisW/libraddask).
 
 
 
 

