#include <cthalassa/cthalassa.hpp>

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

    void Propagator::propagate(const std::vector<double> &times, const std::vector<double> &stateIn, std::vector<std::vector<double>> &statesOut) const {
        // Create copies of the parameter structures
        cthalassa::internal::THALASSAPropagatorStruct settings = settings_;
        cthalassa::internal::THALASSAObjectStruct spacecraft = spacecraft_;

        // State structs
        cthalassa::internal::THALASSAStateStruct rv1, rv2;

        // Determine the number of outputs
        size_t ntime = times.size();

        // Initialise the output vector
        statesOut = std::vector<std::vector<double>>(ntime);

        // Set the initial state
        rv1.mjd = times[0];
        std::copy(stateIn.begin(), stateIn.end(), rv1.RV);
        rv2 = rv1;

        // Add initial state to the output
        statesOut[0] = stateIn;

        // Variable for storing step size
        double tspan;

        // Iterate through the times
        for (size_t itime = 0; itime < times.size() - 1; itime++) {
            // Update times
            rv1.mjd = times[itime];
            tspan = times[itime + 1] - times[itime];
            settings.tspan = tspan;
            settings.tstep = tspan;

            // Update starting state from previous step
            std::copy(std::begin(rv2.RV), std::end(rv2.RV), std::begin(rv1.RV));

            // Propagate step
            cthalassa::internal::thalassa_run(&rv1, &rv2, &spacecraft, &settings);

            // Save state
            statesOut[itime + 1] = std::vector<double>(std::begin(rv2.RV), std::end(rv2.RV));
        }
    };

} // namespace cthalassa