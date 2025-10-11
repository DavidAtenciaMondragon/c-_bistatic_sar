#ifndef RADAR_PARAMS_H
#define RADAR_PARAMS_H

#include <vector>

// Structure for radar transmitter parameters
struct RadarTx {
    double FreqPortadora;      // Carrier frequency [Hz]
    double NumVoltasEsp;       // Number of spiral turns
    double RaioMenorEsp;       // Minor radius of spiral [m]  
    double RaioMaiorEsp;       // Major radius of spiral [m]
    double AltMaiorEsp;        // Maximum altitude of spiral [m]
    double AltMenorEsp;        // Minimum altitude of spiral [m]
    double Vt;                 // Tangential velocity [m/s]
    double PRF;                // Pulse repetition frequency [Hz]
    double lamb;               // Wavelength [m] (calculated)
    double fs;                 // Sampling frequency [Hz]
    double PotenciaTx;         // Transmit power [W]
    double GananciaAnt;        // Antenna gain [dB]
    double FreqMenor;          // Lower frequency [Hz]
    double FreqMaior;          // Upper frequency [Hz]
    double DuracaoPulso;       // Pulse duration [s]
    double AperturaAzimut;     // Azimuth aperture [deg]
    double AperturaElev;       // Elevation aperture [deg]
    double Yaw;                // Yaw angle [deg]
    double InclElev;           // Elevation inclination [deg]
    double Temperatura;        // Temperature [K]
};

// Structure for radar receiver parameters  
struct RadarRx {
    double NumVoltasEsp;       // Number of spiral turns
    double RaioMenorEsp;       // Minor radius of spiral [m]
    double RaioMaiorEsp;       // Major radius of spiral [m]
    double AltMaiorEsp;        // Maximum altitude of spiral [m]
    double AltMenorEsp;        // Minimum altitude of spiral [m]
    double Vt;                 // Tangential velocity [m/s]
    double PRF;                // Pulse repetition frequency [Hz]
    double fs;                 // Sampling frequency [Hz]
    double GananciaAnt;        // Antenna gain [dB]
    double AperturaAzimut;     // Azimuth aperture [deg]
    double AperturaElev;       // Elevation aperture [deg]
    double Yaw;                // Yaw angle [deg]
    double InclElev;           // Elevation inclination [deg]
    double Temperatura;        // Temperature [K]
    double B;                  // Bandwidth [Hz]
};

// Structure for target parameters
struct Target {
    double pos[3];                 // Target position [x, y, z] in meters
    double rcs;                    // Radar cross section
};

// Structure for system parameters
struct System {
    double VelocidadeLuz;          // Speed of light [m/s]
    int IndiceMaximo;              // Maximum index for data matrix
    int NumeroPulsos;              // Number of pulses
    bool PlotarGrafico;            // Whether to plot graphics
    double RefGeografico[3];       // Geographic reference [x, y, z]
    double xMin, xMax;             // X range [m]
    double yMin, yMax;             // Y range [m]
    double discGrid;               // Grid discretization
    double ratioUp;                // Upsampling ratio
};

// Structure for environment parameters
struct Environment {
    double n1;                     // Refractive index 1 (usually air = 1)
    double n2;                     // Refractive index 2 (usually ground = 4)
};

// Structure for processed DEM data
struct ProcessedDEM {
    std::vector<double> X_vec;     // X vector (first row)
    std::vector<double> Y_vec;     // Y vector (first column) 
    std::vector<std::vector<double>> Z_DEM; // Z elevation matrix
};

// Structure for grid processing parameters
struct GridToProc {
    std::vector<double> xAxis;     // X axis coordinates [m]
    std::vector<double> yAxis;     // Y axis coordinates [m] 
    std::vector<double> zAxis;     // Z axis coordinates [m]
};

#endif // RADAR_PARAMS_H