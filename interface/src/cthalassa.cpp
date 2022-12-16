#include <cthalassa/cthalassa.hpp>

#include <math.h>

#include <iostream>

namespace cthalassa {

    Propagator::Propagator(const PropagatorModel &model, const PropagatorPaths &paths, const PropagatorSettings &settings,
                           const PropagatorSpacecraft &spacecraft)
        : model_(model), paths_(paths), settings_(settings), spacecraft_(spacecraft) {
        // Create temporary structures
        cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
        cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;

        // Open THALASSA interface
        cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
    }

    Propagator::~Propagator() {
        // Close the CTHALASSA interface
        cthalassa::internal::thalassa_close();
    }

    void Propagator::propagate(const double &tStart, const double &tEnd, const double &tStep, const std::vector<double> &stateIn, std::vector<double> &timesOut,
                               std::vector<std::vector<double>> &statesOut) const {
        // Create copies of the parameter structures
        cthalassa::internal::THALASSAPropagatorStruct settings = settings_;
        cthalassa::internal::THALASSAObjectStruct spacecraft = spacecraft_;

        // Calculate time span, and number of output times
        double tSpan = tEnd - tStart;
        size_t ntime = std::ceil(tSpan / tStep) + 1;

        // Initialise output vectors
        timesOut = std::vector<double>(ntime);
        statesOut = std::vector<std::vector<double>>(ntime, std::vector<double>(6));

        // Set initial state array
        double initialState[6];
        std::copy(stateIn.begin(), stateIn.end(), initialState);

        // Set propagator times
        settings.tspan = tSpan;
        settings.tstep = tStep;

        // Declare output pointer
        double *output;

        // Propagate
        cthalassa::internal::thalassa_run(&tStart, &initialState[0], &output, &spacecraft, &settings);

        // Extract output
        for (size_t itime = 0; itime < ntime; itime++) {
            // Extract times
            timesOut[itime] = output[itime];

            // Extract states
            for (size_t istate = 0; istate < 6; istate++) {
                statesOut[itime][istate] = output[ntime + itime + istate * ntime];
            }
        }

        // Free output pointer
        free(output);
    }

} // namespace cthalassa