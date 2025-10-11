#include "radar_init.h"

// Function to initialize radar transmitter parameters
// In a real implementation, these would be loaded from JSON files
RadarTx initializeRadarTx() {
    RadarTx strRadarTx;
    
    // Real parameters from radarTx_espiral.json
    strRadarTx.FreqPortadora = 10e9;       // 10 GHz
    strRadarTx.NumVoltasEsp = 2;           // Number of spiral turns
    strRadarTx.RaioMenorEsp = 172.5;       // Minor radius [m]
    strRadarTx.RaioMaiorEsp = 147.5;       // Major radius [m]
    strRadarTx.AltMaiorEsp = 120.0;        // Maximum altitude [m]
    strRadarTx.AltMenorEsp = 80.0;         // Minimum altitude [m]
    strRadarTx.Vt = 120.0;                 // Tangential velocity [m/s]
    strRadarTx.PRF = 400.0;                // PRF [Hz]
    strRadarTx.fs = 60e6;                  // Sampling frequency [Hz]
    strRadarTx.PotenciaTx = 100e3;         // Transmit power [W]
    strRadarTx.GananciaAnt = 6;            // Antenna gain [dB]
    strRadarTx.FreqMenor = -25e6;          // Lower frequency [Hz]
    strRadarTx.FreqMaior = 25e6;           // Upper frequency [Hz] 
    strRadarTx.DuracaoPulso = 1.3e-6;      // Pulse duration [s]
    strRadarTx.AperturaAzimut = 14;        // Azimuth aperture [deg]
    strRadarTx.AperturaElev = 50;          // Elevation aperture [deg]
    strRadarTx.Yaw = 0;                    // Yaw angle [deg]
    strRadarTx.InclElev = 20;              // Elevation inclination [deg]
    strRadarTx.Temperatura = 300;          // Temperature [K]
    
    return strRadarTx;
}

// Function to initialize radar receiver parameters  
RadarRx initializeRadarRx() {
    RadarRx strRadarRx;
    
    // Real parameters from radarRx_espiral.json
    strRadarRx.NumVoltasEsp = 2;           // Number of spiral turns
    strRadarRx.RaioMenorEsp = 172.5;       // Minor radius [m]
    strRadarRx.RaioMaiorEsp = 147.5;       // Major radius [m]
    strRadarRx.AltMaiorEsp = 120.0;        // Maximum altitude [m]
    strRadarRx.AltMenorEsp = 80.0;         // Minimum altitude [m]
    strRadarRx.Vt = 120.0;                 // Tangential velocity [m/s]
    strRadarRx.PRF = 400.0;                // PRF [Hz]
    strRadarRx.fs = 60e6;                  // Sampling frequency [Hz]
    strRadarRx.GananciaAnt = 6;            // Antenna gain [dB]
    strRadarRx.AperturaAzimut = 14;        // Azimuth aperture [deg]
    strRadarRx.AperturaElev = 50;          // Elevation aperture [deg]
    strRadarRx.Yaw = 180;                  // Yaw angle [deg]
    strRadarRx.InclElev = 20;              // Elevation inclination [deg]
    strRadarRx.Temperatura = 300;          // Temperature [K]
    strRadarRx.B = 60e6;                   // Bandwidth [Hz]
    
    return strRadarRx;
}

// Function to initialize target parameters
Target initializeTarget() {
    Target strTarget;
    
    // Initialize target at origin (0,0,0) - position will be set later using updatePointToProc
    strTarget.pos[0] = 0.0;     // X position [m] - initialized at origin
    strTarget.pos[1] = 0.0;     // Y position [m] - initialized at origin
    strTarget.pos[2] = 0.0;     // Z position [m] - initialized at origin
    strTarget.rcs = 1.0;        // Radar cross section
    
    return strTarget;
}

// Function to initialize system parameters
System initializeSystem() {
    System strSystem;
    
    // Real system parameters from system_espiral.json
    strSystem.VelocidadeLuz = 3e8;           // Speed of light [m/s]
    strSystem.IndiceMaximo = 256;            // Maximum index for data matrix
    strSystem.NumeroPulsos = 10000;          // Number of pulses
    strSystem.PlotarGrafico = true;          // Whether to plot graphics
    strSystem.RefGeografico[0] = 0;          // Geographic reference X
    strSystem.RefGeografico[1] = 0;          // Geographic reference Y
    strSystem.RefGeografico[2] = 0;          // Geographic reference Z
    strSystem.xMin = -150;                   // X minimum [m]
    strSystem.xMax = 150;                    // X maximum [m]
    strSystem.yMin = -150;                   // Y minimum [m] 
    strSystem.yMax = 150;                    // Y maximum [m]
    strSystem.discGrid = 1;                  // Grid discretization
    strSystem.ratioUp = 2;                   // Upsampling ratio
    
    return strSystem;
}

// Function to update target position
void updatePointToProc(Target& target, double x, double y, double z) {
    target.pos[0] = x;
    target.pos[1] = y;
    target.pos[2] = z;
}