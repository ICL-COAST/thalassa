#ifndef CTHALASSA
#define CTHALASSA

#include <stddef.h>

/**
 * @brief Struct for physical model parameters
 * 
 * @author Max Hallgarten La Casta
 * @date 2022-11-29
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
 * @date 2022-11-29
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
 * @date 2022-11-29
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
 * @date 2022-11-29
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
 * @date 2022-11-29
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
 * @date 2022-11-29
 *
 * @param[in] model Physical model parameters
 * @param[in] paths Physical model filepaths
 */
void thalassa_open(THALASSAPhysicalModelStruct* model, THALASSAPathStruct* paths);

/**
 * @brief Close the THALASSA interface
 * 
 * @author Max Hallgarten La Casta
 * @date 2022-11-29
 */
void thalassa_close();

/**
 * @brief Propagate an object using THALASSA
 * 
 * @author Max Hallgarten La Casta
 * @date 2022-11-29
 * 
 * @param[in] initialstate Initial state of the object
 * @param[out] finalstate Final state of the object
 * @param[in] object Object parameters
 * @param[in] propagator Propagator parameters
 */
void thalassa_run(THALASSAStateStruct* initialstate, THALASSAStateStruct* finalstate, THALASSAObjectStruct* object, THALASSAPropagatorStruct* propagator);

#endif