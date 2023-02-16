#include <cthalassa/cthalassa.hpp>

#include <cmath>

#ifdef CTHALASSA_USE_FORK
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace cthalassa {

    // Set default static values for PropagatorInstances
    size_t PropagatorInstances::instances_ = 0;
    std::mutex PropagatorInstances::instancesMutex_;
    std::mutex PropagatorInstances::propagationMutex_;
    std::shared_mutex PropagatorInstances::propagationSharedMutex_;

    Propagator::Propagator(const Model &model, const Paths &paths, const Settings &settings, const Spacecraft &spacecraft)
        : model_(model), paths_(paths), settings_(settings), spacecraft_(spacecraft) {
        // Take ownership of the instances mutex
        std::scoped_lock<std::mutex> lock(instancesMutex_);

        // Open the THALASSA interface if only one instance of Propagator exists
        if (instances_ == 1) {
            // Take unique ownership of the shared propagation mutex to ensure that there are no pending propagations
            std::unique_lock<std::shared_mutex> lockShared(propagationSharedMutex_);

            // Make local copies of the model and path structures
            cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
            cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;

            // Open the THALASSA interface
            cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
        }
    }

    Propagator::~Propagator() {
        // Take ownership of the instances mutex
        std::scoped_lock<std::mutex> lock(instancesMutex_);

        // Close the THALASSA interface if only one instance of Propagator remains
        if (instances_ == 1) {
            // Take unique ownership of the shared propagation mutex to ensure that there are no pending propagations
            /// @warning The destructor will block the main thread if there are pending propagations
            std::unique_lock<std::shared_mutex> lockShared(propagationSharedMutex_);

            // Close the THALASSA interface
            cthalassa::internal::thalassa_close();
        }
    }

    void Propagator::propagate(const double &tStart, const double &tEnd, const double &tStep, const std::vector<double> &stateIn, std::vector<double> &timesOut,
                               std::vector<std::vector<double>> &statesOut) const {
        // Take shared ownership of the shared propagation mutex to prevent the THALASSA interface from being closed before all propagations are complete
        std::shared_lock<std::shared_mutex> lockShared(propagationSharedMutex_);

#ifndef CTHALASSA_USE_FORK
        // Take ownership of the propagation mutex
        std::scoped_lock<std::mutex> lock(propagationMutex_);
#endif

        // Create copies of the parameter structures
        cthalassa::internal::THALASSAPropagatorStruct settings = settings_;
        cthalassa::internal::THALASSAObjectStruct spacecraft = spacecraft_;

        // Calculate time span, and number of output times
        double tSpan = tEnd - tStart;
        size_t ntime = std::ceil(tSpan / tStep) + 1;

        // Initialise output vectors
        timesOut.resize(ntime);
        statesOut.resize(ntime, std::vector<double>(6));

        // Set initial state array
        double initialState[6];
        std::copy(stateIn.begin(), stateIn.end(), initialState);

        // Set propagator times
        settings.tspan = tSpan;
        settings.tstep = tStep;

#ifdef CTHALASSA_USE_FORK
        // Declare shared output pointer
        size_t sharedMemorySize = 7 * ntime * sizeof(double);
        double *output = (double *)mmap(NULL, sharedMemorySize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);

        // Propagate
        if (fork() == 0) { // child process
            // Declare local output pointer
            double *outputLocal;

            // Propagate
            cthalassa::internal::thalassa_run(&tStart, &initialState[0], &outputLocal, &spacecraft, &settings);

            // Copy local output to shared output
            std::memcpy(output, outputLocal, sharedMemorySize);

            // Free local output pointer
            free(outputLocal);

            // Exit child process
            exit(0);
        } else { // parent process
            // Wait for child process to finish
            wait(NULL);
        }
#else
        // Declare output pointer
        double *output;

        // Propagate
        cthalassa::internal::thalassa_run(&tStart, &initialState[0], &output, &spacecraft, &settings);
#endif

        // Extract output
        for (size_t itime = 0; itime < ntime; itime++) {
            // Extract times
            timesOut[itime] = output[itime];

            // Extract states
            for (size_t istate = 0; istate < 6; istate++) {
                statesOut[itime][istate] = output[ntime + itime + istate * ntime];
            }
        }

#ifdef CTHALASSA_USE_FORK
        // Free output pointer
        munmap(output, sharedMemorySize);
#endif
    }

} // namespace cthalassa