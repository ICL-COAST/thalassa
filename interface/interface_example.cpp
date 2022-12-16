#include <string.h>

#include <cthalassa/cthalassa.hpp>

int main() {
    cthalassa::PropagatorModel model;

    cthalassa::PropagatorPaths paths = {
        "./data/physical_constants.txt",
        "./data/earth_potential/GRIM5-S1.txt",
        "./data/kernels_to_load.furnsh"
    };

    cthalassa::PropagatorSettings settings;

    cthalassa::PropagatorSpacecraft spacecraft = {
        +8500.000000000000E+00, +13.00000000000000E+00,
        +13.00000000000000E+00, +2.200000000000000E+00,
        +1.500000000000000E+00
    };

    std::vector<double> times = {59043.00000000000E+00, 5.904400000000000E+04, 59043.00000000000E+00 + 365.25};
    std::vector<double> stateIn = {7.000000000000000E+03, 0.000000000000000E+00, 0.000000000000000E+00, 0.000000000000000E+00, 7.546053286792848E+00, 0.000000000000000E+00};
    std::vector<std::vector<double>> statesOut;

    cthalassa::Propagator propagator(model, paths, settings, spacecraft);

    propagator.propagate(times, stateIn, statesOut);

    return 0;
}