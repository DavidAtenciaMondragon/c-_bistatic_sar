#ifndef RADAR_INIT_H
#define RADAR_INIT_H

#include "radar_params.h"

// Function declarations for radar initialization
RadarTx initializeRadarTx();
RadarRx initializeRadarRx();
Target initializeTarget();
System initializeSystem();

// Function to update target position
void updatePointToProc(Target& target, double x, double y, double z);

#endif // RADAR_INIT_H