#ifndef MTHALASSA_PARSERS_H_
#define MTHALASSA_PARSERS_H_

#include "mex.h"

#include <cthalassa/cthalassa.h>

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

/**
 * @brief Parse state struct from Matlab
 *
 * @author Max Hallgarten La Casta
 *
 * @param[in] stateArray Input state struct
 * @param[out] state Output state struct
 */
void parse_state(const mxArray *stateArray, THALASSAStateStruct *state);

#endif