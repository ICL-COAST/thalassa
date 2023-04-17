% Add path for MEX library
addpath("../../lib/")

% Set initial state
state.mjd = 59043.0;
state.RV = [7000, 0, 0, 0, 8, 0];

% Set model parameters
parameters.model.insgrav = 1;
parameters.model.isun    = 1;
parameters.model.imoon   = 1;
parameters.model.idrag   = 1;
parameters.model.iF107   = 1;
parameters.model.iSRP    = 2;
parameters.model.iephem  = 2;
parameters.model.gdeg    = 2;
parameters.model.gord    = 2;

% Set filepaths
parameters.paths.phys_path   = '../../data/physical_constants.txt';
parameters.paths.earth_path  = '../../data/earth_potential/GRIM5-S1.txt';
parameters.paths.kernel_path = '../../data/kernels_to_load.furnsh';

% Set propagator settings
parameters.settings.tol    = 1e-8;
parameters.settings.tspan  = 365.25;
parameters.settings.tstep  = 1.0;
parameters.settings.mxstep = 1e6;
parameters.settings.imcoll = 1;
parameters.settings.eqs    = 2;

% Set spacecraft
parameters.spacecraft.mass      = 1500;
parameters.spacecraft.area_drag = 15;
parameters.spacecraft.area_srp  = 15;
parameters.spacecraft.cd        = 2.2;
parameters.spacecraft.cr        = 1.5;

% Execute THALASSA via the MEX interface
[times, states] = mthalassa(state, parameters);