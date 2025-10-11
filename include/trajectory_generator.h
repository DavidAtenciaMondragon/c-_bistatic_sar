#ifndef TRAJECTORY_GENERATOR_H
#define TRAJECTORY_GENERATOR_H

#include "radar_params.h"
#include "funcao_espiral.h"

/**
 * @brief Generate transmitter trajectory using spiral function
 * @param tx Transmitter radar parameters
 * @return TrajectoryPoints for transmitter
 */
TrajectoryPoints generateTransmitterTrajectory(const RadarTx& tx);

/**
 * @brief Generate receiver trajectory using spiral function with coordinate transformation
 * @param rx Receiver radar parameters
 * @return TrajectoryPoints for receiver (with PxR = -PxR; PyR = -PyR transformation applied)
 */
TrajectoryPoints generateReceiverTrajectory(const RadarRx& rx);

#endif // TRAJECTORY_GENERATOR_H