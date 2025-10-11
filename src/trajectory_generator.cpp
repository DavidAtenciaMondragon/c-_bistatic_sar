#include "trajectory_generator.h"
#include <cstddef>  // For size_t

TrajectoryPoints generateTransmitterTrajectory(const RadarTx& tx) {
    // Generate transmitter trajectory using spiral function
    TrajectoryPoints trajTx = funcao_espiral(tx.NumVoltasEsp, 
                                           tx.RaioMenorEsp,    // Minor radius passed as 2nd param 
                                           tx.RaioMaiorEsp,    // Major radius passed as 3rd param
                                           tx.AltMaiorEsp, 
                                           tx.AltMenorEsp, 
                                           tx.Vt, 
                                           tx.PRF);
    
    return trajTx;
}

TrajectoryPoints generateReceiverTrajectory(const RadarRx& rx) {
    // Generate receiver trajectory using spiral function
    TrajectoryPoints trajRx = funcao_espiral(rx.NumVoltasEsp, 
                                           rx.RaioMenorEsp,    // Minor radius passed as 2nd param
                                           rx.RaioMaiorEsp,    // Major radius passed as 3rd param
                                           rx.AltMaiorEsp, 
                                           rx.AltMenorEsp, 
                                           rx.Vt, 
                                           rx.PRF);

    // Apply transformation for receiver (PxR = -PxR; PyR = -PyR;)
    for (size_t i = 0; i < trajRx.Px.size(); i++) {
        trajRx.Px[i] = -trajRx.Px[i];
        trajRx.Py[i] = -trajRx.Py[i];
    }
    
    return trajRx;
}