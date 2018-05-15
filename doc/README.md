# Introduction to THALASSA
THALASSA is a Fortran orbit propagation code for bodies in the Earth-Moon-Sun system. It works by integrating either Newtonian equations in Cartesian coordinates or regularized equations of motion with the LSODAR (Livermore Solver for Ordinary Differential equations with Automatic Rootfinding).

THALASSA is a command-line tool, and has been developed and tested on MacOS and Ubuntu Linux platforms, using the ``gfortran`` compiler. Attention has been paid to avoid using advanced Fortran constructs: while they greatly improve the capabilities of the language, their portability is usually awful.

The repository also includes some Python3 scripts to perform batch propagations. This feature is currently experimental, but it shouldn't be too difficult for a Python user to generalize the scripts to perform batch propagations on discrete grids of orbital elements.

# THALASSA User Guide
THALASSA reads settings and initial conditions from two text files. Their paths can be specified as arguments to the THALASSA executable,

    ./thalassa.x [settings_file] [object_file]

If no arguments are provided, THALASSA will read settings and initial conditions from the files `in/input.txt` and `in/object.txt`. Sample files are included in the repository.

Fortran is quite strict when interpreting text. When changing the input files, make sure to respect *all* columns and *do not* insert additional rows without changing the code. Explanatory comments in the files are preceded by pound characters.

## Initial conditions file
The initial conditions file contains the initial orbital elements and physical characteristics of the object to be propagated. The orbital elements are expressed in the EMEJ2000 frame; all the epochs are in Terrestrial Time [(Montenbruck and Gill, 2000)](#Montenbruck2000).

In addition to the initial orbital elements, the user also assigns the object mass, areas used in the calculation of atmospheric drag and SRP accelerations, coefficient of drag, and radiation pressure coefficient.

## Settings file
The settings file is divided in four sections.

### Physical model
The first section contains integer flags that allow the user to tune the parameters of the physical model used for the propagation. The meaning of the flags and their allowed values are:
*  `insgrav`: 1 toggles non-spherical gravity, 0 otherwise
*  `isun`: 1 toggles gravitational perturbation from the Sun, 0 otherwise
*  `imoon`: 1 toggles gravitational perturbation from the Moon, 0 otherwise
*  `idrag`: 1 toggles atmospheric drag perturbation, 0 otherwise
*  `iSRP`: 1 toggles solar radiation pressure perturbation, 0 otherwise
*  `iephem`: select the source of ephemerides of the Sun and the Moon. 1 uses SPICE-read ephemerides (DE431 by default), 2 uses a set of simplified analytical ephemerides by [Meeus (1998)](#Meeus1998).
*  `gdeg`: selects the degree of the Earth's gravitational field (up to 95, although a maximum degree of 15 is recommended).
*  `gord`: selects the order of the Earth's gravitational field (has to be less than min(`gdeg`,95)).

### Integration
The second section tunes the parameters of the numerical solver, LSODAR [(Radhakrishnan and Hindmarsh, 1993)](#Radakrishnan1993).
*  `tol`: local truncation error tolerance. Value should be between 1E-4 and 1E-15. This is the main parameter affecting the accuracy of integration.
*  `tspan`: time span of the integration, in days.
*  `tstep`: output time step, in days. Note that the time step *does not* affect the integration accuracy.
*  `mxstep`: maximum number of output steps. Users should not change this value apart from exceptional situations.

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

## THALASSA output
THALASSA is launched by executing `thalassa.x` as specified above.

![Launching THALASSA](/uploads/f2f23ecd72642545bd1774f31ca36602/thalassa_instructions.gif)

The code will write the files `cart.dat` and `orbels.dat` to the directory specified by the user. These contain the numerically integrated trajectory in cartesian coordinates and orbital elements respectively, in the EMEJ2000 reference frame.

## References
1.  <a name="Montenbruck2000"></a>Montenbruck, O. and Gill, E. "Satellite Orbits. Models, Methods, and Applications". Springer-Verlag Berlin Heidelberg, 2000.
2.  <a name="Meeus1998"></a>Meeus, J. "Astronomical Algorithms", 2nd Ed. Willmann-Bell, 1998.
3.  <a name="Bau2015"></a>Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E. "Non-singular orbital elements for special perturbations in the two-body problem". Monthly Notices of the Royal Astronomical Society **454**, pp. 2890-2908, 2015.
3.  <a name="Radhakrishnan1993"></a> Radhakrishnan, K. and Hindmarsh, A. C. "Description and use of LSODE, the Livermore Solver for Ordinary Differential Equations". NASA Reference Publication 1327, Lawrence Livermore National Laboratory Report UCRL-ID-113855, 1993.
4.  <a name="Stiefel1971"></a> Stiefel E. G. and Scheifele G. "Linear and Regular Celestial Mechanics". Springer-Verlag New York Heidelberg Berlin, 1971.