# meanie3D                                                                          

This project provides a very generic implementation of mean-shift clustering on NetCDF data sets. 
The data sets need to follow the "CF-Metadata":http://cfconventions.org. The mean-shift algorithm 
is a robust unsupervised clustering algorithm with a wide range of applications. Meanie3D provides 
a multivariate implementation with no limits as to the number of variables that can be used in 
constructing the Featurespace. 

This software supports 2D as well as 3D data sets and comes with many configuration options. The code 
was written as part of my efforts within "OASE HErZ group":http://www.herz-tb1.uni-bonn.de and contains 
some modules that are specific to this application as a result. If interest is there I might remove all 
specific components at one point. Meanie3D is implemented largely as a C++ template library with only a 
few object files (mostly from the numerical recipes library) for performance reason and uses OpenMP for 
parallelisation. 

The software package comes with a single command in the form of a python script 'meanie3D' which makes 
it very easy to use. The detection, tracking and post-processing steps (such as statistics, visualisation 
of tracks, source data and clusters as well as the track graph (a simple, graphical representation of 
a tracking dictionary for obtaining detailed view on the history of individual clusters and tracks) are 
all controlled through a single configuration file. 

## License

Meanie3D comes under [[MIT license]]. 

## Version

The latest stable version is v1.6.1. Versions are tagged.  

## Build Instructions

### Dependencies

Meanie3D comes with a number of dependencies that need to be installed prior to attempting installation:
* Boost 1.56 or better (http://www.boost.org)
* FLANN 1.8.0 or better (http://www.cs.ubc.ca/research/flann/)
* Blitz++ (http://sourceforge.net/projects/blitz/)
* shapelib v1.3+ (http://shapelib.maptools.org)
* NetCDF 4.2 or better (http://www.unidata.ucar.edu/software/netcdf/) *including* the netcdf4-C++ API
* HDF5 (http://www.hdfgroup.org/HDF5/)
* Python 2.5 or better (https://www.python.org)
* numpy (try running @pip install numpy@ or download and install from http://www.numpy.org)
* NetCDF4-python (try running @pip install netCDF4@ or download and install from http://unidata.github.io/netcdf4-python)
* libradolan (http://meteo-ubonn.github.io/radolan/)

The following libraries may be used, if they are switched on with the appropriate flags
* Visualization:
  * Imagemagick (https://imagemagick.org)
  * VTK 7.0 or better (http://www.vtk.org)
  * Visit 3.0.0 or better (https://wci.llnl.gov/simulation/computer-codes/visit)
* Generated source code documentation:
  * doxygen (https://www.doxygen.nl/download.html)

### Compiler Prerequisites

Meanie3D uses OpenMP by default (-DWITH_OPENMP=1) and requires an OpenMP enabled compiler, such as:

* GNU 9.0 or better (Linux, Mac)
* OpenMP/LLVM (Mac) - a clang implementation supporting OpenMP (http://clang-omp.github.io). Install this 
compiler for optimal results on OSX. You can disable the OpenMP implementation adding the 
flags -DWITH_OPENMP=0 to your cmake call.

### Build instructions

Meanie3D uses CMAKE to generate makefiles. You can use CMAKE's abilities to generate IDE files if you prefer. Start 
by cloning the master branch (for an up-to date but possibly unstable version) or one of the stable releases.
```
  git clone https://github.com/meteo-ubonn/meanie3D.git
```
TODO: revise handling of map data
If you want to download the OASE topology and mapdata file for visualisation, you can obtain this file by adding:
```
  git clone http://git.meteo.uni-bonn.de/git/oase-mapdata
```
Create a build directory:
```
  mkdir meanie3D-make
  cd meanie3D-make
  cmake ../meanie3D
```

### Build presets
A number of presets are provided to make the process easier. Those are selected via the -DPRESET=<preset name> 
flag to cmake. (Example: `cmake -DPRESET=dev-all ../meanie3D`). The available presets are:
* `docker` 
  * Code optimizations for build type 'Release'
  * Core functions and python frontend
* `dev-core`
  * Code optimizations for build type 'Debug'
  * Core functions and python frontend
* `dev-vtk`
  * Code optimizations for build type 'Debug'
  * Core functions and python frontend
  * Visualization code
  * Tests 
  * Documentation
* `dev-all`
  * Code optimizations for build type 'Debug'
  * Core functions and python frontend
  * Visualization code
  * Tests 
  * Documentation
  * RADOLAN, Satellite and KONRAD utilities 
* `prod-core`
  * Code optimizations for build type 'Release'
  * Core functions and python frontend.
* `prod-vtk`
  * Code optimizations for build type 'Release'
  * Core functions and python frontend
  * Visualization code
  * Tests 
  * Documentation
* `prod-all`
  * Code optimizations for build type 'Release'
  * Core functions and python frontend
  * Visualization code
  * Tests 
  * Documentation
  * RADOLAN, Satellite and KONRAD utilities 

The term "core functions" refers to the detection, tracking and track evaluation code. For local 
development in most cases, the preset `dev-vtk` is sufficient. For more aggressive optimisations, 
use `prod-vtk`. If you just need the core functions and none of the visuals, choose the `core` sets 
(`dev-core` or `prod-core`). 

### Available build types
In order to switch optimizations on, tell cmake to use the release build type:
```
  cmake -DCMAKE_BUILD_TYPE=Release ../meanie3D
```
Note that there have been some problems on Linux with aggressive optimization and NetCDF. Your mileage may vary. 
You should try the release build in any event, since it speeds up performance a lot. If you observe unexpected 
problems in reading/writing NetCDF files, you may have fallen victim to the problem and revert to standard build 
(leave the -DCMAKE_BUILD_TYPE=Release).  Once all dependencies are successfully resolved, install the product by 
calling the following: 
```
  make install
```

### Options
If you would like to select your own options, you can leave the PRESET parameter and switch things
on and off using the following options:

#### -DWITH_VTK=ON/OFF
Because of the large footprint of the VTK package, the visualization code is disabled by default. While 
visualization is not necessary to run the algorithm, it can be useful to develop your parameters to 
have visual queues as to what is happening. Setting this flag will result in the following changes:
* The `meanie3D-cfm2vtk` binary will be compiled. This tool can visualize netCDF files in Visit/VTK
* In several places in the code, visualizable output for intermediary steps becomes available:
  * Cluster boundaries/outlines
  * Cluster 'modes' (algorithmic center of a cluster)
  * Cluster geometrical centers
  * Weight function used in detection
  * Mean-shift vectors 
  * Visualization of search window for mean shift.
  * Cluster weight function response
  * Individual variables of the original netCDF (by cluster)
  * Cluster tracks

*Important Notes*: Visualization is switched off in the Docker version. The visualiation code uses
libradolan. If you do switch this on, you will be required to install libradolan as well. 

#### -DWITH_OPENMP=ON/OFF
In order to speed the process up, meanie3D uses OpenMP to parallelize it's computation. This option
is switched on by default.

#### -D WITH_TESTS=ON/OFF
Meanie3D has a number of regression tests, that cover the core algorithms and collection classes. This will 
become important to you if you should decide to work on the core algorithms yourself. The tests are a good 
method of making sure you haven't broken anything critical. The unit tests can then be run by calling
```  
  make test
```

#### -DWITH_PYTHON=ON/OFF
In addition to the C++ binaries, the installation also uses pip to put a python package in place. 
This adds an executable simply called `meanie3D`, which is a front-end to the core functions. It allows 
you to put your clustering and tracking parameters down in the form of a configuration file. The entire 
pipeline is handled based on this configuration file. *This is the recommended way to run the software*. 
For details on the configuration file format, see HOWTO.md. 

#### -DWITH_DOCS=ON/OFF
In order to start developing your own Meanie3D code, it might be useful to have API documentation of the 
various classes in the project. If you have doxygen installed, you can call the following make command to 
create a browsable HTML documentation in doc/html (open the file index.html). 
```
  make docs
```

#### -DWITH_RADOLAN_UTILS=ON/OFF
This will result in compilation of the `meanie3D-radolan2cfm` utility, which converts files in RADOLAN 
format to a cf-metadata compliant netCDF file, which then can be used to run the tracking. 

#### -DWITH_SATELLITE_UTILS=ON/OFF
The package comes with binaries to perform some conversion on satellite data. Those binaries were 
provided in the context of research work for the OASE project. The following binaries will be provided 
if this flag is set:
* `meanie3D-satconv` - Converts spectral radiance to equivalent brightness temperature or vice versa
* `meanie3D-parallax_correction` - Applies parallax correction to mseviri satellite data in OASE composite files.

#### -DWITH_KONRAD_UTILS=ON/OFF
The package comes with a tool `meanie3D-trackstats-conrad` which analyses KONRAD tracks in a way that makes
the data comparable to meanie3D data. This is a specialized tool developed in the context of the OASE project.

## Frequently Asked Quesions 

### I'm having trouble compiling the code
The code has been tested on MacOS X and Suse Linux. Feel free to report problems using the 
"New Issue":http://git.meteo.uni-bonn.de/projects/meanie3d/issues/new link above. Please 
prepend "[BUILD PROBLEM]:" to your subject line. 

### NetCDF/HDF5 trouble
Versions - make sure the NetCDF version is 4.2 or better. Also make certain that you have a HDF5 
installation to match it. If you want to make sure, uninstall HDF5 and NetCDF and re-install them 
from source (HDF5 first, then NetCDF) Commonly there is a problem with the NetCDF C++ API - NetCDF 
has a legacy C++ API from NetCDF3, which does not work with Meanie3D. Make sure you have the correct 
one installed (http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-cxx4-4.2.tar.gz) 

### CMAKE does not find your installed packages
If cmake has trouble locating packages you installed, check the files in the cmake_modules directory 
and adjust the search paths there instead of hardwiring paths into the CMakeLists.txt. If your system 
layout is not considered and you have to make adjustments to the search paths, please send a copy of 
your Find<XZY>.cmake file to juergen.simon@uni-bonn.de so I can include your changes in the next release 
and maybe save someone else the trouble.

## Visualisation
The C++ binaries produce certain files when VTK is enabled. In order to generate imagery or tracking
movies, meanie3D relies on executing python scripts in VisIt (). You will have to install this software
to get access to the visualization. For details on what is possible by scripting VisIt with python,
see the following tutorial: https://www.visitusers.org/index.php?title=VisIt-tutorial-Python-scripting

## More information
Check out the Main Wiki Page at http://git.meteo.uni-bonn.de/projects/meanie3d/wiki for details on how the 
method works, file formats and what the individual binaries do. 