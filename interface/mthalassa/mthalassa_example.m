% Add path for MEX library
addpath("../../lib/")

% Set initial state
state.mjd = 59043.0;
state.RV = [7000, 0, 0, 0, 8, 0];

% Set model parameters
model.insgrav = 1;
model.isun    = 1;
model.imoon   = 1;
model.idrag   = 1;
model.iF107   = 1;
model.iSRP    = 2;
model.iephem  = 2;
model.gdeg    = 2;
model.gord    = 2;

% Set filepaths
paths.phys_path   = '../../data/physical_constants.txt';
paths.earth_path  = '../../data/earth_potential/GRIM5-S1.txt';
paths.kernel_path = '../../data/kernels_to_load.furnsh';

% Set propagator settings
settings.tol    = 1e-8;
settings.tspan  = 365.25;
settings.tstep  = 1.0;
settings.mxstep = 1e6;
settings.imcoll = 1;
settings.eqs    = 2;

% Set spacecraft
spacecraft.mass      = 1500;
spacecraft.area_drag = 15;
spacecraft.area_srp  = 15;
spacecraft.cd        = 2.2;
spacecraft.cr        = 1.5;

% Execute THALASSA via the MEX interface
[times, states] = mthalassa(state, model, paths, settings, spacecraft);