#ifndef CTHALASSA
#define CTHALASSA

#include <stddef.h>

/**
 * @brief Struct for physical model parameters
 *
 * @author Max Hallgarten La Casta
 */
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

/**
 * @brief Struct for propagator parameters
 *
 * @author Max Hallgarten La Casta
 */
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

/**
 * @brief Struct for object parameters
 *
 * @author Max Hallgarten La Casta
 */
typedef struct THALASSAObjectStruct {
    // Physical characteristics
    double mass;
    double area_drag;
    double area_srp;
    double cd;
    double cr;
} THALASSAObjectStruct;

/**
 * @brief Struct for state vectors
 *
 * @author Max Hallgarten La Casta
 */
typedef struct THALASSAStateStruct {
    // Epoch
    double mjd;

    // State
    double RV[6];
} THALASSAStateStruct;

/**
 * @brief Struct for model filepaths
 *
 * @author Max Hallgarten La Casta
 */
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

/**
 * @brief Open the THALASSA interface
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] model Physical model parameters
 * @param[in] paths Physical model filepaths
 */
void thalassa_open(THALASSAPhysicalModelStruct *model, THALASSAPathStruct *paths);

/**
 * @brief Close the THALASSA interface
 *
 * @author Max Hallgarten La Casta
 */
void thalassa_close();

/**
 * @brief Execute a propagation with THALASSA
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] initialtime Initial times
 * @param[in] initialstate Initial Cartesian states
 * @param[out] outputmatrix Matrix containing propagation results
 * @param[out] object Object parameters
 * @param[out] propagator Propagator parameters
 */
void thalassa_run(const double *initialtime, const double *initialstate, double **outputmatrix, const THALASSAObjectStruct *object,
                  const THALASSAPropagatorStruct *propagator);

#endif