#include <math.h>
#include <string.h>

#include "mex.h"

#include <cthalassa/cthalassa.h>

void parse_model(const mxArray *modelArray, THALASSAPhysicalModelStruct *model) {
    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(modelArray);

    // Check for the correct number of fields
    if (nfields != 9) {
        mexErrMsgTxt("Incorrect number of fields in physical model structure");
    }

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(modelArray, ifield);

        // Update fields
        if (strcmp(fname, "insgrav") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->insgrav = round(tmp);
        } else if (strcmp(fname, "isun") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->isun = round(tmp);
        } else if (strcmp(fname, "imoon") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->imoon = round(tmp);
        } else if (strcmp(fname, "idrag") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->idrag = round(tmp);
        } else if (strcmp(fname, "iF107") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iF107 = round(tmp);
        } else if (strcmp(fname, "iSRP") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iSRP = round(tmp);
        } else if (strcmp(fname, "iephem") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iephem = round(tmp);
        } else if (strcmp(fname, "gdeg") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->gdeg = round(tmp);
        } else if (strcmp(fname, "gord") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->gord = round(tmp);
        } else {
            mexErrMsgTxt("Unknown parameter in physical model structure");
        }
    }
}

void parse_paths(const mxArray *pathsArray, THALASSAPathStruct *paths) {
    // Pointers for parsing
    char *tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(pathsArray);

    // Check for the correct number of fields
    if (nfields != 3) {
        mexErrMsgTxt("Incorrect number of fields in path structure");
    }

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(pathsArray, ifield);

        // Update fields
        if (strcmp(fname, "phys_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->phys_path, tmp);
            paths->phys_path_len = strlen(tmp);
        } else if (strcmp(fname, "earth_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->earth_path, tmp);
            paths->earth_path_len = strlen(tmp);
        } else if (strcmp(fname, "kernel_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->kernel_path, tmp);
            paths->kernel_path_len = strlen(tmp);
        } else {
            mexErrMsgTxt("Unknown parameter in paths structure");
        }
    }
}

void parse_propagator(const mxArray *settingsArray, THALASSAPropagatorStruct *settings) {
    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(settingsArray);

    // Check for the correct number of fields
    if (nfields != 6) {
        mexErrMsgTxt("Incorrect number of fields in settings structure");
    }

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(settingsArray, ifield);

        // Update fields
        if (strcmp(fname, "tol") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tol = tmp;
        } else if (strcmp(fname, "tspan") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tspan = tmp;
        } else if (strcmp(fname, "tstep") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tstep = tmp;
        } else if (strcmp(fname, "mxstep") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->mxstep = tmp;
        } else if (strcmp(fname, "imcoll") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->imcoll = round(tmp);
        } else if (strcmp(fname, "eqs") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->eqs = round(tmp);
        } else {
            mexErrMsgTxt("Unknown parameter in settings structure");
        }
    }
}

void parse_spacecraft(const mxArray *spacecraftArray, THALASSAObjectStruct *spacecraft) {
    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(spacecraftArray);

    // Check for the correct number of fields
    if (nfields != 5) {
        mexErrMsgTxt("Incorrect number of fields in spacecraft structure");
    }

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(spacecraftArray, ifield);

        // Update fields
        if (strcmp(fname, "mass") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->mass = tmp;
        } else if (strcmp(fname, "area_drag") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->area_drag = tmp;
        } else if (strcmp(fname, "area_srp") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->area_srp = tmp;
        } else if (strcmp(fname, "cd") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->cd = tmp;
        } else if (strcmp(fname, "cr") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->cr = tmp;
        } else {
            mexErrMsgTxt("Unknown parameter in spacecraft structure");
        }
    }
}

void parse_state(const mxArray *stateArray, THALASSAStateStruct *state) {
    // Pointers for parsing
    double tmp;
    double *tmpArray;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(stateArray);

    // Check for the correct number of fields
    if (nfields != 2) {
        mexErrMsgTxt("Incorrect number of fields in state structure");
    }

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(stateArray, ifield);

        // Update fields
        if (strcmp(fname, "mjd") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(stateArray, 0, ifield));
            state->mjd = tmp;
        } else if (strcmp(fname, "RV") == 0) {
            tmpArray = mxGetPr(mxGetFieldByNumber(stateArray, 0, ifield));
            memcpy(state->RV, tmpArray, sizeof(state->RV));
        } else {
            mexErrMsgTxt("Unknown parameter in state structure");
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Extract inputs
    const mxArray *stateArray = prhs[0];
    const mxArray *modelArray = prhs[1];
    const mxArray *pathsArray = prhs[2];
    const mxArray *settingsArray = prhs[3];
    const mxArray *spacecraftArray = prhs[4];

    // Declare parameter structures
    THALASSAStateStruct state;
    THALASSAPhysicalModelStruct model;
    THALASSAPathStruct paths;
    THALASSAPropagatorStruct settings;
    THALASSAObjectStruct spacecraft;

    // Parse parameter structures
    parse_state(stateArray, &state);
    parse_model(modelArray, &model);
    parse_paths(pathsArray, &paths);
    parse_propagator(settingsArray, &settings);
    parse_spacecraft(spacecraftArray, &spacecraft);

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

    // Open THALASSA interface
    thalassa_open(&model, &paths);

    // Declare pointers for the initial time and state
    const double *inputState = &state.RV[0];

    // Calculate memory size
    size_t memorySize = 6 * ntime * sizeof(double);

    // Allocate memory for output state
    double *outputState = (double *)malloc(memorySize);

    // Propagate
    thalassa_run(&ntime, timesOut, inputState, outputState, &spacecraft, &settings);

    // Copy output
    memcpy(statesOut, outputState, memorySize);

    // Free output matrix
    free(outputState);

    // Close THALASSA interface
    thalassa_close();
}