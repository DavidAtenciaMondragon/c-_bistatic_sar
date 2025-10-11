#ifndef FUNCAO_ESPIRAL_H
#define FUNCAO_ESPIRAL_H

#include <vector>

// Structure to hold the trajectory points
struct TrajectoryPoints {
    std::vector<double> Px;  // X coordinates [m]
    std::vector<double> Py;  // Y coordinates [m] 
    std::vector<double> Pz;  // Z coordinates [m]
};

/**
 * @brief Calculate spiral trajectory points
 * 
 * This function returns the x, y and z points of a spiral trajectory
 * given the input variables:
 * 
 * @param N Number of desired turns
 * @param R_minor Minor radius of the spiral [m] (passed as 2nd parameter in MATLAB call)
 * @param R_major Major radius of the spiral [m] (passed as 3rd parameter in MATLAB call)
 * @param H_high Maximum height of the spiral [m]
 * @param H_low Minimum height of the spiral [m]
 * @param Vt Tangential velocity [m/s]
 * @param PRF Pulse repetition frequency [Hz]
 * 
 * @return TrajectoryPoints structure containing Px, Py, Pz coordinates
 */
TrajectoryPoints funcao_espiral(double N, double R_minor, double R_major, 
                               double H_high, double H_low, double Vt, double PRF);

#endif // FUNCAO_ESPIRAL_H