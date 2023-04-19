#include <math.h>
#include <string.h>

#include "mex.h"

#include <cthalassa/cthalassa.h>

#include <mthalassa/mthalassa.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Extract inputs
    const mxArray *timeArray = prhs[0];
    const mxArray *stateArray = prhs[1];
    const mxArray *parameterArray = prhs[2];

    // Extract times
    size_t ntime;
    double *times;
    parse_times(timeArray, &ntime, &times);
    
    // Extract initial state
    double *inputState;
    parse_state(stateArray, &inputState);

    // Parse inputs
    THALASSAPhysicalModelStruct model;
    THALASSAPathStruct paths;
    THALASSAPropagatorStruct settings;
    THALASSAObjectStruct spacecraft;
    parse_parameters(parameterArray, &model, &paths, &settings, &spacecraft);

    // Create output matrix
    plhs[0] = mxCreateDoubleMatrix(6, (mwSize)ntime, 0);
    double *statesOut = mxGetDoubles(plhs[0]);

    // Open THALASSA interface
    thalassa_open(&model, &paths);

    // Propagate
    thalassa_run(&ntime, times, inputState, statesOut, &spacecraft, &settings);

    // Close THALASSA interface
    thalassa_close();
}