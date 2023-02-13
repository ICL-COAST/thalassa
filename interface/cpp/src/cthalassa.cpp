#include <cthalassa/cthalassa.hpp>

#include <cmath>

#ifdef CTHALASSA_USE_FORK
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace cthalassa {

    Propagator::Propagator(const PropagatorModel &model, const PropagatorPaths &paths, const PropagatorSettings &settings,
                           const PropagatorSpacecraft &spacecraft, const bool &INTERFACE_ISOLATION)
        : model_(model), paths_(paths), settings_(settings), spacecraft_(spacecraft), INTERFACE_ISOLATION_(INTERFACE_ISOLATION) {
#ifdef CTHALASSA_USE_FORK
        // Open THALASSA interface if interface isolation is not requested
        if (!INTERFACE_ISOLATION_) {
            cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
            cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;
            cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
        }
#else
        // Always open THALASSA interface when not using fork
        cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
        cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;
        cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
#endif
    }

    Propagator::~Propagator() {
#ifdef CTHALASSA_USE_FORK
        // Close THALASSA interface if interface isolation is not requested
        if (!INTERFACE_ISOLATION_) {
            cthalassa::internal::thalassa_close();
        }
#else
        // Always close THALASSA interface when not using fork
        cthalassa::internal::thalassa_close();
#endif
    }

    void Propagator::propagate(const double &tStart, const double &tEnd, const double &tStep, const std::vector<double> &stateIn, std::vector<double> &timesOut,
                               std::vector<std::vector<double>> &statesOut) const {
        // Create copies of the parameter structures
        cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
        cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;
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

#ifdef CTHALASSA_USE_FORK
        // Declare shared output pointer
        size_t sharedMemorySize = 7 * ntime * sizeof(double);
        double *output = (double *)mmap(NULL, sharedMemorySize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);

        // Propagate
        if (fork() == 0) { // child process
            // Open THALASSA interface if isolation is enabled
            if (INTERFACE_ISOLATION_) {
                cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
            }

            // Declare local output pointer
            double *outputLocal;

            // Propagate
            cthalassa::internal::thalassa_run(&tStart, &initialState[0], &outputLocal, &spacecraft, &settings);

            // Copy local output to shared output
            std::memcpy(output, outputLocal, sharedMemorySize);

            // Free local output pointer
            free(outputLocal);

            // Close the THALASSA interface if isolation is enabled
            if (INTERFACE_ISOLATION_) {
                cthalassa::internal::thalassa_close();
            }

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
#else
        // Free output pointer
        free(output);
#endif
    }

} // namespace cthalassa