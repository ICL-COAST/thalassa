#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <cthalassa.h>

THALASSAPathStruct *create_thalassa_path_struct(char *phys_path, char *earth_path, char *kernel_path) {
    // Allocate pointer
    THALASSAPathStruct *paths = malloc(sizeof(THALASSAPathStruct));

    // Physical constants path
    strcpy(paths->phys_path, phys_path);
    paths->phys_path_len = strlen(phys_path);

    // Earth model path
    strcpy(paths->earth_path, earth_path);
    paths->earth_path_len = strlen(earth_path);

    // SPICE kernel path
    strcpy(paths->kernel_path, kernel_path);
    paths->kernel_path_len = strlen(kernel_path);

    // Return pointer
    return paths;
}

int main () {
    // Declare physical model parameters, and paths
    THALASSAPhysicalModelStruct model = {1, 1, 1, 1, 1, 2, 1, 2, 2};
    THALASSAPathStruct *paths = create_thalassa_path_struct("./data/physical_constants.txt", "./data/earth_potential/GRIM5-S1.txt", "./data/kernels_to_load.furnsh");

    // Declare initial and final state
    THALASSAStateStruct initialstate = {
        59043.0,
        {7.000000000000000E+03,
         0.000000000000000E+00,
         0.000000000000000E+00,
         0.000000000000000E+00,
         7.546053286792848E+00,
         0.000000000000000E+00}
    };
    THALASSAStateStruct finalstate;

    // Declare object parameters
    THALASSAObjectStruct object = {8500.0, 13.0, 13.0, 2.2, 1.5};

    // Declare propagator parameters
    THALASSAPropagatorStruct propagator = {1e-13, 1.0, 1.0, 1e6, 0, 2};

    // Initialise THALASSA
    thalassa_open(&model, paths);

    // Run THALASSA
    thalassa_run(&initialstate, &finalstate, &object, &propagator);

    // Close THALASSA
    thalassa_close();

    // Free allocated memory
    free(paths);
}