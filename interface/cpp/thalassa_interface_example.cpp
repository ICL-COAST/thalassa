#include <cthalassa/cthalassa.hpp>

int main() {
    cthalassa::PropagatorModel model;

    cthalassa::PropagatorPaths paths = {"./data/physical_constants.txt", "./data/earth_potential/GRIM5-S1.txt", "./data/kernels_to_load.furnsh"};

    cthalassa::PropagatorSettings settings;

    cthalassa::PropagatorSpacecraft spacecraft = {+8500.000000000000E+00, +13.00000000000000E+00, +13.00000000000000E+00, +2.200000000000000E+00,
                                                  +1.500000000000000E+00};

    cthalassa::Propagator propagator(model, paths, settings, spacecraft);

    std::vector<double> stateIn = {7.000000000000000E+03, 0.000000000000000E+00, 0.000000000000000E+00,
                                   0.000000000000000E+00, 7.546053286792848E+00, 0.000000000000000E+00};

    double tStart = 59043;
    double tEnd = tStart + 7;
    double tStep = 1;

    std::vector<double> timesOut;
    std::vector<std::vector<double>> statesOut;

    propagator.propagate(tStart, tEnd, tStep, stateIn, timesOut, statesOut);

    return 0;
}