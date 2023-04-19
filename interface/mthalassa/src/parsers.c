#include <mthalassa/parsers.h>

#include <math.h>
#include <string.h>

void parse_times(const mxArray *timeArray, size_t *ntime, double **times) {
    // Extract number of time array dimensions
    const mwSize ntimedim = mxGetNumberOfDimensions(timeArray);

    // Ensure that the time array is 2D
    if (ntimedim != 2) {
        mexErrMsgTxt("Incompatible time vector dimensions");
    }

    // Extract time array dimensions
    const mwSize *timedim = mxGetDimensions(timeArray);

    // Extract number of times
    if (timedim[0] == 1) {
        *ntime = timedim[1];
    } else if (timedim[1] == 1) {
        *ntime = timedim[0];
    } else {
        mexErrMsgTxt("Incompatible time vector dimensions");
    }

    // Ensure that there are at least two times
    if (*ntime < 2) {
        mexErrMsgTxt("Insufficient number of times");
    }

    // Update times vector pointer
    *times = mxGetDoubles(timeArray);
}

void parse_state(const mxArray *stateArray, double **state) {
    // Extract number of state array dimensions
    const mwSize nstatedim = mxGetNumberOfDimensions(stateArray);

    // Ensure that the state array is 2D
    if (nstatedim != 2) {
        mexErrMsgTxt("Incompatible state vector dimensions");
    }

    // Extract state array dimensions
    const mwSize *statedim = mxGetDimensions(stateArray);

    // Declare variable for the number of state variables
    size_t nstate;

    // Extract number of state variables
    if (statedim[0] == 1) {
        nstate = statedim[1];
    } else if (statedim[1] == 1) {
        nstate = statedim[0];
    } else {
        mexErrMsgTxt("Incompatible initial state dimensions");
    }

    // Ensure that there are six state variables
    if (nstate != 6) {
        mexErrMsgTxt("Incorrect number of state variables");
    }

    // Update state vector pointer
    *state = mxGetDoubles(stateArray);
}

void parse_parameters(const mxArray *parameterArray, THALASSAPhysicalModelStruct *model, THALASSAPathStruct *paths, THALASSAPropagatorStruct *settings,
                      THALASSAObjectStruct *spacecraft) {
    // Declare the exit code
    int exitcode = 4;

    // Pointers for parsing
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(parameterArray);

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(parameterArray, ifield);

        // Parse structs
        if (strcmp(fname, "model") == 0) {
            parse_model(mxGetFieldByNumber(parameterArray, 0, ifield), model);
            --exitcode;
        } else if (strcmp(fname, "paths") == 0) {
            parse_paths(mxGetFieldByNumber(parameterArray, 0, ifield), paths);
            --exitcode;
        } else if (strcmp(fname, "settings") == 0) {
            parse_propagator(mxGetFieldByNumber(parameterArray, 0, ifield), settings);
            --exitcode;
        } else if (strcmp(fname, "spacecraft") == 0) {
            parse_spacecraft(mxGetFieldByNumber(parameterArray, 0, ifield), spacecraft);
            --exitcode;
        }
    }

    // Throw parsing error
    if (exitcode != 0) {
        mexErrMsgTxt("Incomplete parameters");
    }
}

void parse_model(const mxArray *modelArray, THALASSAPhysicalModelStruct *model) {
    // Declare exit code with number of expected fields
    int exitcode = 9;

    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(modelArray);

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(modelArray, ifield);

        // Update fields
        if (strcmp(fname, "insgrav") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->insgrav = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "isun") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->isun = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "imoon") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->imoon = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "idrag") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->idrag = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "iF107") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iF107 = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "iSRP") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iSRP = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "iephem") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->iephem = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "gdeg") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->gdeg = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "gord") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(modelArray, 0, ifield));
            model->gord = round(tmp);
            --exitcode;
        }
    }

    // Throw parsing error
    if (exitcode != 0) {
        mexErrMsgTxt("Incomplete model");
    }
}

void parse_paths(const mxArray *pathsArray, THALASSAPathStruct *paths) {
    // Declare exit code with number of expected fields
    int exitcode = 3;

    // Pointers for parsing
    char *tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(pathsArray);

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(pathsArray, ifield);

        // Update fields
        if (strcmp(fname, "phys_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->phys_path, tmp);
            paths->phys_path_len = strlen(tmp);
            --exitcode;
        } else if (strcmp(fname, "earth_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->earth_path, tmp);
            paths->earth_path_len = strlen(tmp);
            --exitcode;
        } else if (strcmp(fname, "kernel_path") == 0) {
            tmp = mxArrayToString(mxGetFieldByNumber(pathsArray, 0, ifield));
            strcpy(paths->kernel_path, tmp);
            paths->kernel_path_len = strlen(tmp);
            --exitcode;
        }
    }

    // Throw parsing error
    if (exitcode != 0) {
        mexErrMsgTxt("Incomplete paths");
    }
}

void parse_propagator(const mxArray *settingsArray, THALASSAPropagatorStruct *settings) {
    // Declare exit code with number of expected fields
    int exitcode = 6;

    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(settingsArray);

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(settingsArray, ifield);

        // Update fields
        if (strcmp(fname, "tol") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tol = tmp;
            --exitcode;
        } else if (strcmp(fname, "tspan") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tspan = tmp;
            --exitcode;
        } else if (strcmp(fname, "tstep") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->tstep = tmp;
            --exitcode;
        } else if (strcmp(fname, "mxstep") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->mxstep = tmp;
            --exitcode;
        } else if (strcmp(fname, "imcoll") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->imcoll = round(tmp);
            --exitcode;
        } else if (strcmp(fname, "eqs") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(settingsArray, 0, ifield));
            settings->eqs = round(tmp);
            --exitcode;
        }
    }

    // Throw parsing error
    if (exitcode != 0) {
        mexErrMsgTxt("Incomplete settings");
    }
}

void parse_spacecraft(const mxArray *spacecraftArray, THALASSAObjectStruct *spacecraft) {
    // Declare exit code with number of expected fields
    int exitcode = 5;

    // Pointers for parsing
    double tmp;
    const char *fname;

    // Find number of fields
    const mwSize nfields = mxGetNumberOfFields(spacecraftArray);

    // Iterate through fields
    for (mwSize ifield = 0; ifield < nfields; ifield++) {
        // Extract field name
        fname = mxGetFieldNameByNumber(spacecraftArray, ifield);

        // Update fields
        if (strcmp(fname, "mass") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->mass = tmp;
            --exitcode;
        } else if (strcmp(fname, "area_drag") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->area_drag = tmp;
            --exitcode;
        } else if (strcmp(fname, "area_srp") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->area_srp = tmp;
            --exitcode;
        } else if (strcmp(fname, "cd") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->cd = tmp;
            --exitcode;
        } else if (strcmp(fname, "cr") == 0) {
            tmp = mxGetScalar(mxGetFieldByNumber(spacecraftArray, 0, ifield));
            spacecraft->cr = tmp;
            --exitcode;
        }
    }

    // Throw parsing error
    if (exitcode != 0) {
        mexErrMsgTxt("Incomplete spacecraft");
    }
}
