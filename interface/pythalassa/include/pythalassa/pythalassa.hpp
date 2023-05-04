#include <Eigen/Core>

#include <cthalassa/cthalassa.hpp>

namespace pythalassa {

    /**
     * @brief Wrapper class for PyTHALASSA to use CTHALASSA
     *
     * @author Max Hallgarten La Casta
     */
    class Propagator : public cthalassa::Propagator {

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
        Propagator(const cthalassa::Model &model, const cthalassa::Paths &paths, const cthalassa::Settings &settings, const cthalassa::Spacecraft &spacecraft)
            : cthalassa::Propagator(model, paths, settings, spacecraft) {}

        /**
         * @brief Propagate a state
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] times Times
         * @param[in] stateIn Initial state
         * @return Eigen::MatrixXd Output states
         */
        Eigen::MatrixXd propagate(const Eigen::VectorXd &times, const Eigen::VectorXd &stateIn) const {
            // Declare output
            Eigen::MatrixXd statesOut;

            // Propagate orbit
            cthalassa::Propagator::propagate(times, stateIn, statesOut);

            // Return states
            return statesOut;
        }
    };

} // namespace pythalassa