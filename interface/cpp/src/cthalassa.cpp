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
    Model PropagatorInstances::model_;
    Paths PropagatorInstances::paths_;

    Propagator::Propagator(const Model &model, const Paths &paths, const Settings &settings, const Spacecraft &spacecraft)
        : settings_(settings), spacecraft_(spacecraft) {
        // Take ownership of the instances mutex
        std::unique_lock<std::mutex> lock_instancesMutex(instancesMutex_);

        // Open the THALASSA interface if only one instance of Propagator exists
        if (instances_ == 1) {
            // Unique lock to make changes to the THALASSA interface
            std::unique_lock<std::shared_mutex> lock_propagationSharedMutex(propagationSharedMutex_);

            // Set model and paths
            model_ = model;
            paths_ = paths;

            // Make local copies of the model and path structures
            cthalassa::internal::THALASSAPhysicalModelStruct modelTemp = model_;
            cthalassa::internal::THALASSAPathStruct pathsTemp = paths_;

            // Open the THALASSA interface
            cthalassa::internal::thalassa_open(&modelTemp, &pathsTemp);
        } else {
            // Throw errors if the model or paths do match what is already loaded
            if ((model != model_)) {
                throw std::runtime_error("Requested model does not match the existing model");
            }
            if ((paths != paths_)) {
                throw std::runtime_error("Requested paths do not match the existing paths");
            }
        }
    }

    Propagator::~Propagator() {
        // Unique lock to enable changes to instances
        std::unique_lock<std::mutex> lock_instancesMutex(instancesMutex_);

        // Close the THALASSA interface if only one instance of Propagator remains
        if (instances_ == 1) {
            // Unique lock to make changes to the THALASSA interface
            /// @warning The destructor will block the main thread if there are pending propagations
            std::unique_lock<std::shared_mutex> lock_propagationSharedMutex(propagationSharedMutex_);

            // Close the THALASSA interface
            cthalassa::internal::thalassa_close();
        }
    }

    void Propagator::propagate(const double &tStart, const double &tEnd, const double &tStep, const std::vector<double> &stateIn, std::vector<double> &timesOut,
                               std::vector<std::vector<double>> &statesOut) const {
        // Calculate number of steps
        const size_t Ntimes = std::ceil((tEnd - tStart) / tStep) + 1;

        // Resize time vector
        timesOut.resize(Ntimes);

        // Declare temporary time variable
        double tTemp;

        // Populate time vector
        for (size_t itime = 0; itime < Ntimes; itime++) {
            // Calculate time
            tTemp = tStart + itime * tStep;

            // Clamp number below end time
            if (tTemp > tEnd) {
                tTemp = tEnd;
            }

            // Update time vector
            timesOut[itime] = tTemp;
        }

        // Propagate
        propagate(timesOut, stateIn, statesOut);
    };

    void Propagator::propagate(const std::vector<double> &times, const std::vector<double> &stateIn, std::vector<std::vector<double>> &statesOut) const {
        // Declare shared lock to prevent changes to the THALASSA interface while propagations remain to be executed
        std::shared_lock<std::shared_mutex> lock_propagationSharedMutex{propagationSharedMutex_, std::defer_lock};

        // Declare shared lock to prevent changes to the local settings and spacecraft while propagations remain to be executed
        std::shared_lock<std::shared_mutex> lock_propagationLocalSharedMutex{propagationLocalSharedMutex_, std::defer_lock};

        // Simultaneous shared lock
        std::lock(lock_propagationSharedMutex, lock_propagationLocalSharedMutex);

#ifndef CTHALASSA_USE_FORK
        // Unique lock to ensure serial execution
        std::unique_lock<std::mutex> lock_propagationMutex(propagationMutex_);
#endif

        // Create copies of the parameter structures
        cthalassa::internal::THALASSAPropagatorStruct settings = settings_;
        cthalassa::internal::THALASSAObjectStruct spacecraft = spacecraft_;

        // Calculate time span
        double tStart = times.front();
        double tEnd = times.back();
        double tSpan = tEnd - tStart;

        // Extract number of output states
        const int ntime = times.size();

        // Calculate mean time step
        double tStep = (tEnd - tStart) / (ntime - 1);

        // Initialise output vectors
        statesOut.resize(ntime, std::vector<double>(6));

#ifdef CTHALASSA_USE_FORK
        // Declare shared output pointer
        size_t sharedMemorySize = 6 * ntime * sizeof(double);
        double *output = (double *)mmap(NULL, sharedMemorySize, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANON, -1, 0);

        // Fork process
        pid_t pid = fork();

        // Branch based on process identifier
        if (pid == 0) { // child process
            // Declare local output pointer
            double *outputLocal;

            // Propagate
            cthalassa::internal::thalassa_run(&ntime, &times.front(), &stateIn.front(), &spacecraft, &settings, &outputLocal);

            // Copy local output to shared output
            std::memcpy(output, outputLocal, sharedMemorySize);

            // Free local output pointer
            free(outputLocal);

            // Exit child process
            exit(0);
        } else { // parent process
            // Wait for the child process to finish
            waitpid(pid, NULL, 0);
        }
#else
        // Declare output pointer
        double *output;

        // Propagate
        cthalassa::internal::thalassa_run(&ntime, &times.front(), &stateIn.front(), &spacecraft, &settings, &output);
#endif

        // Extract output
        for (size_t itime = 0; itime < ntime; itime++) {
            // Extract states
            for (size_t istate = 0; istate < 6; istate++) {
                statesOut[itime][istate] = output[itime + istate * ntime];
            }
        }

#ifdef CTHALASSA_USE_FORK
        // Free output pointer
        munmap(output, sharedMemorySize);
#endif
    }

    Model Propagator::getModel() const {
        // Return model
        return model_;
    }

    Paths Propagator::getPaths() const {
        // Return paths
        return paths_;
    }

    void Propagator::setSettings(const Settings &settings) {
        // Take unique ownership of the local shared propagation mutex to ensure that there are no pending propagations
        std::unique_lock<std::shared_mutex> lock(propagationLocalSharedMutex_);

        // Update settings
        settings_ = settings;
    }

    Settings Propagator::getSettings() {
        // Take shared ownership of the local shared propagation mutex to ensure that the settings are not changed while read
        std::shared_lock<std::shared_mutex> lock(propagationLocalSharedMutex_);

        // Return settings
        return settings_;
    }

    void Propagator::setSpacecraft(const Spacecraft &spacecraft) {
        // Take unique ownership of the local shared propagation mutex to ensure that there are no pending propagations
        std::unique_lock<std::shared_mutex> lock(propagationLocalSharedMutex_);

        // Update spacecraft parameters
        spacecraft_ = spacecraft;
    }

    Spacecraft Propagator::getSpacecraft() {
        // Take shared ownership of the local shared propagation mutex to ensure that the spacecraft is not changed while read
        std::shared_lock<std::shared_mutex> lock(propagationLocalSharedMutex_);

        // Return spacecraft
        return spacecraft_;
    }

} // namespace cthalassa