#include <math.h>
#include <string.h>

#include "mex.h"

#include <cthalassa/cthalassa.h>

#include <mthalassa/mthalassa.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Extract inputs
    const mxArray *stateArray = prhs[0];
    const mxArray *parameterArray = prhs[1];

    // Declare parameter structures
    THALASSAStateStruct state;
    THALASSAPhysicalModelStruct model;
    THALASSAPathStruct paths;
    THALASSAPropagatorStruct settings;
    THALASSAObjectStruct spacecraft;

    // Parse inputs
    parse_state(stateArray, &state);
    parse_parameters(parameterArray, &model, &paths, &settings, &spacecraft);

    // Calculate number of times
    const size_t ntime = ceil(settings.tspan / settings.tstep) + 1;

    // Extract start and end times
    double tStart = state.mjd;
    double tEnd = state.mjd + settings.tspan;

    // Declare temporary time variable
    double tTemp;

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)ntime, 0);
    plhs[1] = mxCreateDoubleMatrix(6, (mwSize)ntime, 0);

    // Declare pointers to output matrices
    double *timesOut = mxGetDoubles(plhs[0]);
    double *statesOut = mxGetDoubles(plhs[1]);

    // Populate time vector
    for (size_t itime = 0; itime < ntime; itime++) {
        // Calculate time
        tTemp = tStart + itime * settings.tstep;

        // Clamp number below end time
        if (tTemp > tEnd) {
            tTemp = tEnd;
        }

        // Update time vector
        timesOut[itime] = tTemp;
    }

    // Declare pointers for the initial time and state
    const double *inputState = &state.RV[0];

    // Calculate memory size
    size_t memorySize = 6 * ntime * sizeof(double);

    // Allocate memory for output state
    double *outputState = (double *)malloc(memorySize);

    // Open THALASSA interface
    thalassa_open(&model, &paths);

    // Propagate
    thalassa_run(&ntime, timesOut, inputState, outputState, &spacecraft, &settings);

    // Close THALASSA interface
    thalassa_close();

    // Copy output
    memcpy(statesOut, outputState, memorySize);

    // Free output matrix
    free(outputState);
}