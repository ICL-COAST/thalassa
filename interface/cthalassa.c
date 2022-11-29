#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct THALASSAPhysicalModelStruct {
    // Physical model
    int insgrav;
    int isun;
    int imoon;
    int idrag;
    int iF107;
    int iSRP;
    int iephem;
    int gdeg;
    int gord;
} THALASSAPhysicalModelStruct;

typedef struct THALASSAPropagatorStruct {
    // Integration
    double tol;
    double tspan;
    double tstep;
    double mxstep;
    int imcoll;

    // Equations of motion
    int eqs;
} THALASSAPropagatorStruct;

typedef struct THALASSAObjectStruct {
    // Physical characteristics
    double mass;
    double area_drag;
    double area_srp;
    double cd;
    double cr;
} THALASSAObjectStruct;

typedef struct THALASSAStateStruct {
    // Epoch
    double mjd;

    // State
    double RV[6];
} THALASSAStateStruct;

typedef struct THALASSAPathStruct {
    // Physical constants path
    char phys_path[512];
    size_t phys_path_len;

    // Earth model path
    char earth_path[512];
    size_t earth_path_len;

    // SPICE kernel path
    char kernel_path[512];
    size_t kernel_path_len;  
} THALASSAPathStruct;

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

void thalassa_open(THALASSAPhysicalModelStruct*, THALASSAPathStruct*);
void thalassa_close();
void thalassa_run(THALASSAStateStruct*, THALASSAStateStruct*, THALASSAObjectStruct*, THALASSAPropagatorStruct*);

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