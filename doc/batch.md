# Batch mode
In the `batch` directory you will find Python3 scripts to launch THALASSA over a grid of epochs and initial conditions, which are defined through orbital elements in the EMEJ2000 reference frame.

## Grid creation
To create a grid of orbital elements, edit the file `griddef.json`. The file specifies the extremes of each dimension in the 7-dimensional grid, and the points for each dimension.
Once this is specified, launch `./coegrid.py [outDir]`, where `outDir` is the directory where the output trajectories will be saved.
The script displays the total number of propagations that will be prepared, and asks for confirmation to the user before creating a `grid.dat` file in the output directory, which is just a text table in which the initial conditions for each point of the grid are saved.
Along with `grid.dat`, the script also saves the `griddef.json` to the output directory in order to ensure full reproducibility.

### Directory structure
Each set of initial conditions is assigned a sequential ID (SID) in the grid table. The order in which the initial conditions are saved is the inverse of the usual order of the orbital elements, i.e. $M, \omega, \Omega, i, e, a, \text{MJD}$.
The script also creates in `outDir` the following directory structure:

    +-- outDir
        +-- C001
        |   +-- S0000000001
        |   +-- S0000000002
        |   +-- ...
        |   +-- S0000010000
        |
        +-- C002
        |   +-- S0000010001
        |   +-- ...
        |   +-- S0000020000
        |
        +-- ...
        |
        +-- CYYY
The subdirectories `SXXXXXXXXXX`, where `XXXXXXXXXX` is the 10-digit SID, contain the `input.txt`, `object.txt` and the resulting trajectories for each of the propagation.
To avoid extremely long dirtrees that could case problems on some file systems, they are chunked into the directories `CYYY`.
By default, each `CYYY` directory contains 10000 `SXXXXXXXXXX` subdirectory.
     
**WARNING!!** The number of propagations can quickly get _big_. Think carefully about the size of your grid and the available computational and storage resources before starting a batch.