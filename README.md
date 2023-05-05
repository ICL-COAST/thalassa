# THALASSA
THALASSA is a Fortran orbit propagation code for bodies in the Earth-Moon-Sun system. It works by integrating either Newtonian equations in Cartesian coordinates or regularized equations of motion with the LSODAR (Livermore Solver for Ordinary Differential equations with Automatic Root-finding).

THALASSA is a command-line tool, and has been developed and tested on MacOS and Ubuntu Linux platforms, using the ``gfortran`` compiler. Attention has been paid to avoid using advanced Fortran constructs: while they greatly improve the capabilities of the language, their portability has been found to be problematic. This might change in the future.

The repository also includes some Python3 scripts to perform batch propagations. This feature is currently experimental, but it shouldn't be too difficult for a Python user to generalize the scripts to perform batch propagations on discrete grids of orbital elements. Interfaces for C/C++ (CTHALASSA), Matlab (MTHALASSA), and Python (PyTHALASSA) have been developed to avoid the file input/output bottleneck.

Details on the mathematical fundamentals of THALASSA are contained in [Amato et al., 2018](#Amato2018).

# Giving credit
If you use THALASSA in your research, please consider giving credit by citing the [article](https://doi.org/10.1007/s10569-019-9897-1) detailing the mathematical fundamentals (Amato et al., 2019). In addition, you can cite the THALASSA [ASCL entry](http://ascl.net/1905.018) according to the [ASCL guidelines](http://ascl.net/home/getwp/351).

# Getting Started
## Prerequisites
THALASSA was originally developed using the `gfortran` compiler. Compilation of CTHALASSA, MTHALASSA, and PyTHALASSA require additional compilers from the GNU compiler collection. Additionally, compilation of MTHALASSA requires a Matlab installation with an activated licence.

Compilers:
* `gcc`
* `g++`
* `gfortran`

Build system:
* `CMake`
* `make` (or an alternative generator for CMake)

The source files for the `sofa` and `spice` libraries are automatically downloaded and extracted during CMake's configuration stage. `Eigen` is an optional dependency for the CTHALASSA interface, and a required dependency for the PyTHALASSA interface.

THALASSA is also dependent on `spice` kernels if precise ephemerides are required. A bash script (`./dependencies.sh`) is provided for convenience to download THALASSA's default `spice` kernels.

## Compilation
### Traditional
THALASSA is compiled with two `CMake` commands. These configure the project and then build it:
```
cmake -B ./build -S .
```
```
cmake --build ./build
```
By default, all targets will be built, including THALASSA, CTHALASSA, and MTHALASSA.

#### Known Issues/Limitations

* If `CMake` freezes during the configuration stage, and Matlab is installed, this may be due to the Matlab license not being activated, preventing the successful configuration of the `MEX` compiler.

* THALASSA was designed for single threaded execution. It is *not* thread safe by default. CTHALASSA spawns subprocesses for each propagation with `fork` to provide thread safety, however this depends on POSIX compliance. This can be deactivated, if required, by appending `-DCTHALASSA_USE_FORK=OFF` to the `CMake` configuration command. This is deactivated on non-UNIX systems.

#### Compiling MTHALASSA on Apple Machines with ARM
The current release of Matlab for MacOS uses `x86_64` via the Rosetta translation environment. Consequently, the `MEX` compiler targets `x86_64`. However, by default `CMake` will target the native architecture of the machine, in this case `arm64`. Furthermore, older versions of `CMake` are not aware of Rosetta installations, and will attempt to use the ARM version of the `MEX` compiler.

To compile a compatible version of MTHALASSA, the follow actions must be taken:
* Install an updated version of `CMake` (>=3.26)
* Install `x86_64` versions of the C, C++, and Fortran compilers 

The following flags must be appended to the `CMake` configuration command to force the architecture, and force the use of the `x86_64` compilers:
```
-DCMAKE_OSX_ARCHITECTURES=x86_64
```
```
-DCMAKE_C_COMPILER=<C compiler path>
```
```
-DCMAKE_CXX_COMPILER=<C++ compiler path>
```
```
-DCMAKE_Fortran_COMPILER=<Fortran compiler path>
```
This will build `x86_64` versions of all of the THALASSA components, ensuring that it can be used via Matlab.

An example of installing `x86_64` versions of the GNU compilers via `Homebrew` is provided below:
```
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
```
arch -x86_64 /usr/local/bin/brew install gcc
```
The compilers should be available under `/usr/local/Cellar/gcc/<version>/bin/`.

### Container
THALASSA can take advantage of containerisation by being built with the `Dockerfile` in the repository's root directory:
```
docker build -t thalassa .
```
This will automatically download and install the required compile-time and run-time dependencies. The kernels required by `spice` (listed in `./data/kernels_to_load.furnsh`) must be downloaded manually, and provided to the container via a volume bind mount.

By default, the Docker container will be based on Debian Bullseye, however an Alpine-based version is available by adding `--build-arg BASE_DISTRO=alpine` to the build command.

Note: testing has revealed that the Alpine-based image produces different solutions to Debian-based machines and images. It is suspected that this is due to differences in the dependencies. Users are advised to proceed with caution.

## Usage
### THALASSA
THALASSA reads settings and initial conditions from two text files. Their paths can be specified as arguments to the THALASSA executable:
```
./thalassa_main [settings_file] [object_file]
```
If no arguments are provided, THALASSA will read settings and initial conditions from the files `./in/input.txt` and `./in/object.txt`. Sample files are included in the repository.

Fortran is quite strict when interpreting text. When changing the input files, make sure to respect *all* columns and *do not* insert additional rows without changing the code. Explanatory comments in the files are preceded by pound characters.

The code will write the files `cart.dat` and `orbels.dat` to the directory specified by the user. These contain the numerically integrated trajectory in cartesian coordinates and orbital elements respectively, in the EMEJ2000 reference frame. Additionally, the code writes a file `stats.dat` containing integration statistics along with the value of the orbital elements at the final epoch, and a `propagation.log` file which contains diagnostics for the current propagation.

You should check the `stats.dat` file for any errors that might have taken place during the execution. In particular, THALASSA will write to the log file and to stdout the exit code of the propagation. This is an integer number with the following meaning:
* `0`: nominal exit, reached end time specified in `input.txt`
* `1`: nominal exit, collided with the Earth (atmospheric re-entry).
* `-2`: floating point exception, detected NaNs in the state vector. This is usually caused by the orbit having become hyperbolic when using EDromo, or in some more exotic cases due to a solver tolerance that's too large.
* `-3`: maximum number of steps reached. Try specifying a larger time step in the input file.
* `-10`: unknown exception, try debugging to check what's the problem.

An example of using the container is provided below, including the binds for the input, output, and kernel directories:
```
docker run \
       --rm \
       -v <host input directory>:/thalassa/in/:ro \
       -v <host output directory>:/thalassa/out/ \
       -v <host kernel directory>:/thalassa/data/kernels/:ro \
       thalassa
```

Alternatively, individual files within the required directories can be bound, if one does not want to expose the entire directory to the container.

Note: the output directory specified in `input.txt` MUST be `./out/` to ensure that the output of THALASSA is saved successfully through the bind to the host.

### CTHALASSA
CTHALASSA is available as a `CMake` target for easy integration with other `CMake`-based projects. This automatically includes the header files, and takes care of its dependency on THALASSA.

Alternatively, for more manual use, the library is available in `./lib/`, and the header files in `./interface/cthalassa/include/`.

It is recommended to use the definitions in `cthalassa.hpp` as these are more user friendly, and take care of pre- and post-propagation processing. Nevertheless, the C definitions used by CTHALASSA to interface with Fortran are available in `cthalassa.h`.

CTHALASSA can be used by creating a `cthalassa::Propagator` object:
```
cthalassa::Propagator propagator(model, paths, settings, spacecraft);
```
The constructor takes multiple structures as input, including model parameters, filepaths to the physical model files (as described below), propagator settings, and spacecraft properties. These structures were designed to shadow the original text files used by THALASSA.

The propagator can be called with the `propagate` method:
```
propagator.propagate(tStart, tEnd, tStep, stateIn, timesOut, statesOut);
```
or
```
propagator.propagate(times, stateIn, statesOut);
```

An example for using CTHALASSA is available (`./interface/cthalassa/cthalassa_example.cpp`) which propagates an orbit via the interface.

Note: it is recommended to use absolute filepaths for the physical model files. Furthermore, it is recommended to provide a custom `.furnsh` file also containing absolute filepaths when `spice` ephemerides are being used.

### MTHALASSA
The recommended method for using MTHALASSA is to add THALASSA's library directory to Matlab's path at the beginning of a script:
```
addpath(<THALASSA lib directory>)
```

MTHALASSA can then be used with the following signature:
```
[statesOut] = mthalassa(times, stateIn, parameters)
```

A documentation file is automatically generated which can be accessed from Matlab:
```
help mthalassa
```

An example of using MTHALASSA is available (`./interface/mthalassa/mthalassa_example.m`) which propagates an orbit via the interface.

Note: it is recommended to use absolute filepaths for the physical model files. Furthermore, it is recommended to provide a custom `.furnsh` file also containing absolute filepaths when `spice` ephemerides are being used.

### PyTHALASSA
The recommended method for using PyTHALASSA is to add THALASSA's library directory to Python's path at the beginning of a script:
```
import os
import sys
sys.path.append(<THALASSA lib directory>)
import pythalassa
```

PyTHALASSA can be used by creating a `pythalassa.Propagator` object:
```
propagator = pythalassa.Propagator(model, paths, settings, spacecraft)
```
The constructor takes multiple structures as input, including model parameters, filepaths to the physical model files (as described below), propagator settings, and spacecraft properties. These structures were designed to shadow the original text files used by THALASSA.

The propagator can be called with the `propagate` method:
```
statesOut = propagator.propagate(times, stateIn)
```

An example of using PyTHALASSA is available (`./interface/pythalassa/pythalassa_example.py`) which propagates an orbit via the interface.

Note: it is recommended to use absolute filepaths for the physical model files. Furthermore, it is recommended to provide a custom `.furnsh` file also containing absolute filepaths when `spice` ephemerides are being used.

# THALASSA Input Format
## Initial conditions file
The initial conditions file contains the initial orbital elements and physical characteristics of the object to be propagated. The orbital elements are expressed in the EMEJ2000 frame; all the epochs are in Terrestrial Time [(Montenbruck and Gill, 2000)](#Montenbruck2000).

In addition to the initial orbital elements, the user also assigns the object mass, areas used in the calculation of atmospheric drag and SRP accelerations, coefficient of drag, and radiation pressure coefficient.

## Settings file
The settings file is divided in four sections.

### Physical model
The first section contains integer flags that allow the user to tune the parameters of the physical model used for the propagation. The meaning of the flags and their allowed values are:
*  `insgrav`: 1 toggles non-spherical gravity, 0 otherwise
*  `isun`: values >1 are interpreted as the order of the Legendre expansion for the solar gravitational perturbing acceleration. 1 toggles the acceleration using the full expression in rectangular coordinates. 0 disables the perturbation. See [Amato et al. (2018)](#Amato2018) for details.
*  `imoon`: values >1 are interpreted as the order of the Legendre expansion for the lunar gravitational perturbing acceleration. 1 toggles the acceleration using the full expression in rectangular coordinates. 0 disables the perturbation. See [Amato et al. (2018)](#Amato2018) for details.
*  `idrag`: select atmospheric model. 0 = no drag, 1 = patched exponential model [(table 8-4 of Vallado and McClain, 2013)](#Vallado2013), 2 = US 1976 Standard Atmosphere [(NASA et al., 1976)](#US1976), 3 = Jacchia 1977 [(Jacchia, 1977)](#Jacchia1977), 4 = NRLMSISE-00 [(Picone et al., 2000)](#Picone2000).
*  `iF107`: 1 toggles variable F10.7 flux, 0 uses the constant value specified in `./data/physical_constants.txt`
*  `iSRP`: select SRP model. 0 = no SRP, 1 = SRP with no eclipses, 2 = SRP with conical shadow using the $\nu$ factor from [Montenbruck and Gill (2000)](#Montenbruck2000).
*  `iephem`: select the source of ephemerides of the Sun and the Moon. 1 uses SPICE-read ephemerides (DE431 by default), 2 uses a set of simplified analytical ephemerides by [Meeus (1998)](#Meeus1998).
*  `gdeg`: selects the degree of the Earth's gravitational field (up to 95, although a maximum degree of 15 is recommended).
*  `gord`: selects the order of the Earth's gravitational field (has to be less than min(`gdeg`,95)).

### Integration
The second section tunes the parameters of the numerical solver, LSODAR [(Radhakrishnan and Hindmarsh, 1993)](#Radakrishnan1993).
*  `tol`: local truncation error tolerance. Value should be between 1E-4 and 1E-15. This is the main parameter affecting the accuracy of integration.
*  `tspan`: time span of the integration, in days.
*  `tstep`: output time step, in days. Note that the time step *does not* affect the integration accuracy.
*  `mxstep`: maximum number of output steps. Users should not change this value apart from exceptional situations.
*  `imcol`: 1 toggles the check for collision with the Moon, 0 disables it. Activating the check roughly doubles the CPU time needed to compute lunar perturbations, therefore activate the check only if necessary.

### Equations of motion
The third section only contains the `eqs` flag, which selects the set of equations of motion to be integrated. The value of `eqs` corresponds to the following equations:
1.  Cowell formulation (Newtonian equations in Cartesian coordinates)
2.  EDromo orbital elements, including physical time as a dependent variable [(Baù et al., 2015)](#Bau2015)
3.  EDromo orbital elements, including the constant time element as a dependent variable
4.  EDromo orbital elements, including the linear time element as a dependent variable
5.  Kustaanheimo-Stiefel coordinates, including the physical time as a dependent variable [(section 9 of Stiefel and Scheifele, 1971)](#Stiefel1971)
6.  Kustaanheimo-Stiefel coordinates, including the linear time element as a dependent variable
7.  Stiefel-Scheifel orbital elements, including the physical time as a dependent variable [(section 19 of Stiefel and Scheifele, 1971)](#Stiefel1971)
8.  Stiefel-Scheifel orbital elements, including the linear time element as a dependent variable

The choice of equations depends on the type of orbit being integrated. For LEOs and MEOs, sets 2 to 6 are recommended. Sets 2 to 4 are particularly efficient for HEOs but should be avoided if there's any chance for the orbital energy to change sign, as the integration will fail in such a case.
As a rule of thumb, weakly-perturbed orbits can be integrated most efficiently by using orbital elements.
Strongly-perturbed orbits should be integrated using coordinates, i.e. sets 1, 5, 6.

### Output
The last section contains settings for the output of THALASSA.
*  `verb`: 1 toggles the printing of the current propagation progress, 0 otherwise
*  `out`:  Full path to the directory where THALASSA saves the output files. The path **always** starts at column 5, and ends with a `/`.

It is recommended to untoggle the `verb` flag if THALASSA is used to propagate batches of initial conditions. Failure to do so could unnecessarily clutter `stdout`.

## Physical data files
The directory `./data/` stores files containing information on the physical model used by THALASSA. `./data/earth_potential/GRIM5-S1.txt` contains the harmonic coefficients of the GRIM5-S1 potential, in its native format.
The file `./data/physical_constants.txt` contains several astronomical and physical constants that are used during the integration.
The constant values of the solar flux and of the planetary index and amplitude are specified here, along with the height at which the orbiter is considered to have re-entered the atmosphere of the Earth.

# References
1.  <a name="Amato2018"></a>Amato, D., Bombardelli, C., Baù, G., Morand, V., and Rosengren, A. J. "Non-averaged regularized formulations as an alternative to semi-analytical orbit propagation methods". Submitted to Celestial Mechanics and Dynamical Astronomy, 2018.
2.  <a name="Montenbruck2000"></a>Montenbruck, O., and Gill, E. "Satellite Orbits. Models, Methods, and Applications". Springer-Verlag Berlin Heidelberg, 2000.
3.  <a name="Vallado2013"></a>Vallado, D. A., and McClain, W. D. "Fundamentals of Astrodynamics and Applications". Microcosm Press, 2013.
4. <a name="US1976"></a> NASA, NOAA, and US Air Force, "U.S. Standard Atmosphere, 1976". Technical Report NASA-TM-X-74335, October 1976.
5. <a name="Jacchia1977"></a> Jacchia, L. G. "Thermospheric Temperature, Density, and Composition: New Models". SAO Special Report, **375**, 1977.
6. <a name="Picone2000"></a> Picone, J. M., Hedin, A. E., Drob, D .P., and Aikin, A. C. "NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues". Journal of Geophysical Research: Space Physics, **107**(A12):15–1–16, 2002.
7.  <a name="Meeus1998"></a>Meeus, J. "Astronomical Algorithms", 2nd Ed. Willmann-Bell, 1998.
8.  <a name="Bau2015"></a>Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E. "Non-singular orbital elements for special perturbations in the two-body problem". Monthly Notices of the Royal Astronomical Society **454**, pp. 2890-2908, 2015.
9.  <a name="Radhakrishnan1993"></a> Radhakrishnan, K. and Hindmarsh, A. C. "Description and use of LSODE, the Livermore Solver for Ordinary Differential Equations". NASA Reference Publication 1327, Lawrence Livermore National Laboratory Report UCRL-ID-113855, 1993.
10. <a name="Stiefel1971"></a> Stiefel E. G. and Scheifele G. "Linear and Regular Celestial Mechanics". Springer-Verlag New York Heidelberg Berlin, 1971.