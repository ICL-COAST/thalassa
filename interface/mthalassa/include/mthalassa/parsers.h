#ifndef MTHALASSA_PARSERS_H_
#define MTHALASSA_PARSERS_H_

#include "mex.h"

#include <cthalassa/cthalassa.h>

/**
 * @brief Parse times from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] timeArray Input time array
 * @param[out] ntime Number of times
 * @param[out] times Time vector
 */
void parse_times(const mxArray *timeArray, size_t *ntime, double **times);

/**
 * @brief Parse state from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] stateArray Input state array
 * @param[out] state State vector
 */
void parse_state(const mxArray *stateArray, double **state);

/**
 * @brief Parse parameters struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] parameterArray Input parameters struct
 * @param[out] model Output model struct
 * @param[out] paths Output paths struct
 * @param[out] settings Output settings struct
 * @param[out] spacecraft Output spacecraft struct
 */
void parse_parameters(const mxArray *parameterArray, THALASSAPhysicalModelStruct *model, THALASSAPathStruct *paths, THALASSAPropagatorStruct *settings,
                      THALASSAObjectStruct *spacecraft);

/**
 * @brief Parse model struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] modelArray Input model struct
 * @param[out] model Output model struct
 */
void parse_model(const mxArray *modelArray, THALASSAPhysicalModelStruct *model);

/**
 * @brief Parse paths struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] pathsArray Input paths struct
 * @param[out] paths Output paths struct
 */
void parse_paths(const mxArray *pathsArray, THALASSAPathStruct *paths);

/**
 * @brief Parse propagator struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] settingsArray Input settings struct
 * @param[out] settings Output settings struct
 */
void parse_propagator(const mxArray *settingsArray, THALASSAPropagatorStruct *settings);

/**
 * @brief Parse spacecraft struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] spacecraftArray Input spacecraft struct
 * @param[out] spacecraft Output spacecraft struct
 */
void parse_spacecraft(const mxArray *spacecraftArray, THALASSAObjectStruct *spacecraft);

#endif