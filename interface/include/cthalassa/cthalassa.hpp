#ifndef CTHALASSA_HPP_
#define CTHALASSA_HPP_

#include <cstring>
#include <string>
#include <vector>

namespace cthalassa::internal {
    extern "C" {
    // clang-format off
        #include <cthalassa/cthalassa.h>
    // clang-format on
    }
} // namespace cthalassa::internal

namespace cthalassa {

    /// @brief Gravity models
    enum PropagatorModelGravity { SPHERICAL = 0, NONSPHERICAL = 1 };

    /// @brief Sun model
    enum PropagatorModelSun { SUN_DISABLED = 0, SUN_ENABLED = 1 };

    /// @brief Moon model
    enum PropagatorModelMoon { MOON_DISABLED = 0, MOON_ENABLED = 1 };

    /// @brief Drag model
    enum PropagatorModelDrag { DRAG_DISABLED = 0, DRAG_WERTZ = 1, DRAG_US76 = 2, DRAG_J77 = 3, DRAG_NRLMSISE00 = 4 };

    /// @brief Flux model
    enum PropagatorModelFlux { FLUX_CONSTANT = 0, FLUX_VARIABLE = 1 };

    /// @brief Solar Radiation Pressure (SRP) model
    enum PropagatorModelSRP { SRP_DISABLED = 0, SRP_ENABLED = 1, SRP_ENABLED_CONICAL = 2 };

    /// @brief Ephemerides
    enum PropagatorModelEphemerides { EPHEM_DE431 = 1, EPHEM_SIMPLE = 2 };

    /**
     * @brief Structure containing propagator model settings
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct PropagatorModel {
        /// @brief Gravity model
        PropagatorModelGravity insgrav = NONSPHERICAL;

        /// @brief Sun model
        PropagatorModelSun isun = SUN_ENABLED;

        /// @brief Moon model
        PropagatorModelMoon imoon = MOON_ENABLED;

        /// @brief Drag model
        PropagatorModelDrag idrag = DRAG_WERTZ;

        /// @brief Flux model
        PropagatorModelFlux iF107 = FLUX_VARIABLE;

        /// @brief Solar Radiation Pressure (SRP) model
        PropagatorModelSRP iSRP = SRP_ENABLED_CONICAL;

        /// @brief Ephemerides
        PropagatorModelEphemerides iephem = EPHEM_SIMPLE;

        /// @brief Gravity model degree
        int gdeg = 2;

        /// @brief Gravity model order
        int gord = 2;

        /// @brief Method to implicitly cast the structure to its equivalent in
        /// CTHALASSA
        operator cthalassa::internal::THALASSAPhysicalModelStruct() const {
            return cthalassa::internal::THALASSAPhysicalModelStruct{insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord};
        }
    } PropagatorModel;

    /**
     * @brief Structure containing filepaths for THALASSA
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct PropagatorPaths {
        /// @brief Filepath for THALASSA physical constants
        std::string phys_path;

        /// @brief Filepath for Earth geopotential model
        std::string earth_path;

        /// @brief Filepath for SPICE kernel
        std::string kernel_path;

        /// @brief Method to implicitly cast the structure to its equivalent in
        /// CTHALASSA
        operator cthalassa::internal::THALASSAPathStruct() const {
            // Declare new path struct
            cthalassa::internal::THALASSAPathStruct paths;

            // Set physical constants path
            strcpy(paths.phys_path, phys_path.c_str());
            paths.phys_path_len = phys_path.size();

            // Set Earth model path
            strcpy(paths.earth_path, earth_path.c_str());
            paths.earth_path_len = earth_path.size();

            // Set SPICE kernel path
            strcpy(paths.kernel_path, kernel_path.c_str());
            paths.kernel_path_len = kernel_path.size();

            // Return path struct
            return paths;
        }
    } PropagatorPaths;

    /// @brief Propagator formulations
    enum PropagatorSettingsEquations { COWELL = 1, EDROMO_T = 2, EDROMO_C = 3, EDROMO_L = 4, KS_T = 5, KS_L = 6, STISCHE_T = 7, STISCHE_L = 8 };

    /**
     * @brief Structure containing propagator settings
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct PropagatorSettings {
        /// @brief Integration tolerance
        double tol = 1.0E-13;

        /// @brief Maximum number of integration/output steps
        double mxstep = 1.0E6;

        /// @brief Moon collision check flag
        int imcoll = 0;

        /// @brief Propagator formulation
        PropagatorSettingsEquations eqs = EDROMO_T;

        /// @brief Propagation time span [solar days]
        double tspan;

        /// @brief Step size [solar days]
        double tstep;

        /// @brief Method to implicitly cast the structure to its equivalent in
        /// CTHALASSA
        operator cthalassa::internal::THALASSAPropagatorStruct() const {
            return cthalassa::internal::THALASSAPropagatorStruct{tol, tspan, tstep, mxstep, imcoll, eqs};
        }
    } PropagatorSettings;

    /**
     * @brief Structure containing spacecraft parameters
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct PropagatorSpacecraft {
        /// @brief Spacecraft mass [kg]
        double mass;

        /// @brief Drag area [m^2]
        double area_drag;

        /// @brief Solar Radiation Pressure (SRP) area [m^2]
        double area_srp;

        /// @brief Coefficient of drag [-]
        double cd;

        /// @brief Coefficient of reflectivity [-]
        double cr;

        /// @brief Method to implicitly cast the structure to its equivalent in
        /// CTHALASSA
        operator cthalassa::internal::THALASSAObjectStruct() const { return cthalassa::internal::THALASSAObjectStruct{mass, area_drag, area_srp, cd, cr}; }
    } PropagatorSpacecraft;

    /**
     * @brief THALASSA orbit propagator
     *
     * @author Max Hallgarten La Casta
     */
    class Propagator {

    private:
        /// @brief Model settings
        const PropagatorModel model_;

        /// @brief Filepaths
        const PropagatorPaths paths_;

        /// @brief Propagator settings
        const PropagatorSettings settings_;

        /// @brief Spacecraft
        const PropagatorSpacecraft spacecraft_;

    public:
        /**
         * @brief Construct a new Propagator object
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] model Model settings
         * @param[in] paths Filepaths
         * @param[in] settings Propagator settings
         * @param[in] spacecraft Spacecraft
         */
        Propagator(const PropagatorModel &model, const PropagatorPaths &paths, const PropagatorSettings &settings, const PropagatorSpacecraft &spacecraft);

        /**
         * @brief Destroy the Propagator object
         *
         * @author Max Hallgarten La Casta
         */
        ~Propagator();

        /**
         * @brief Propagate the spacecraft
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] times Vector of times in Modified Julian Date (MJD)
         * @param[in] stateIn Initial state vector
         * @param[out] statesOut State vectors at the requested times
         */
        void propagate(const std::vector<double> &times, const std::vector<double> &stateIn, std::vector<std::vector<double>> &statesOut) const;
    };

} // namespace cthalassa

#endif