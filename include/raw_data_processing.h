#ifndef RAW_DATA_PROCESSING_H
#define RAW_DATA_PROCESSING_H

#include "radar_params.h"
#include "funcao_espiral.h"
#include "load_dem.h"
#include <vector>
#include <complex>

// Structure for trajectory matrices
struct TrajectoryMatrix {
    std::vector<std::vector<double>> Tx; // Transmitter positions [Nx3]
    std::vector<std::vector<double>> Rx; // Receiver positions [Nx3]
    std::vector<double> P;               // Target position [3x1]
};

// Structure for range calculations
struct RangeData {
    std::vector<double> r1Tx;  // Range from Tx to refraction point
    std::vector<double> r2Tx;  // Range from refraction point to target
    std::vector<double> r1Rx;  // Range from Rx to refraction point  
    std::vector<double> r2Rx;  // Range from target to refraction point
};

// Structure for refraction points (ida = outbound path)
struct RefractionPoints {
    std::vector<std::vector<double>> P_refrac_ida;    // Refraction points for Tx->Target path [Nx3]
    std::vector<std::vector<double>> P_refrac_volta;  // Refraction points for Target->Rx path [Nx3]
};

// Structure for reflection points  
struct ReflectionPoints {
    std::vector<std::vector<double>> P_reflection;   // Reflection points [Nx3]
    // Additional reflection data can be added here
};

// Structure for raw data processing
struct RawDataResult {
    std::vector<std::vector<std::complex<double>>> auxData;  // Auxiliary data matrix
    std::vector<std::vector<std::complex<double>>> rawData;  // Raw data matrix
    std::vector<std::vector<std::complex<double>>> rootData; // Compressed data matrix
    std::vector<std::complex<double>> chirp;                 // CHIRP signal
};

/**
 * @brief Prepare trajectory data for processing
 * 
 * Converts trajectory points to matrix format and prepares target position
 * Equivalent to MATLAB:
 * Tx = [PxT.', PyT.', PzT.'];
 * Rx = [PxR.', PyR.', PzR.'];
 * P  = strTarget.pos.';
 */
TrajectoryMatrix prepareTrajectoryData(const TrajectoryPoints& trajTx, 
                                      const TrajectoryPoints& trajRx,
                                      const Target& target);

// Overloaded version for matrix-based data
TrajectoryMatrix prepareTrajectoryData(const std::vector<std::vector<double>>& Tx,
                                      const std::vector<std::vector<double>>& Rx,
                                      const std::vector<double>& target_pos);

/**
 * @brief Prepare DEM data for processing
 * 
 * Extracts vectors and converts DEM data to processing format
 * Equivalent to MATLAB:
 * strDEM.X_vec = X_DEM(1,:);
 * strDEM.Y_vec = Y_DEM(:,1).';
 * strDEM.Z_DEM = double(Z_DEM);
 */
ProcessedDEM prepareDEMData(const ProcessedDEM& demData);

/**
 * @brief Initialize environment parameters
 * 
 * Sets up refractive indices for ray tracing calculations
 * Equivalent to MATLAB:
 * n1 = 1; n2 = 4;
 * strEnvironment.n1 = n1; strEnvironment.n2 = n2;
 */
Environment initializeEnvironment();

/**
 * @brief Calculate ranges using Fermat's principle for refraction
 * 
 * Implements the full calculaSlantRangeFermat function that computes
 * reflection and refraction points using Fermat's principle for accurate
 * ray tracing through two media with different refractive indices.
 * 
 * Equivalent to MATLAB:
 * [strReflexao, strRefraccoes] = calculaSlantRangeFermat(strDEM,Tx,Rx,P,n1,n2,bPlotVerbose);
 * 
 * @param demData Processed DEM data for surface intersection calculations
 * @param trajectories Trajectory matrices containing Tx, Rx, and target positions
 * @param env Environment parameters with refractive indices n1 and n2
 * @param strReflexao Output: reflection points and associated data
 * @param strRefraccoes Output: refraction points for ida and volta paths
 * @param bPlotVerbose Enable verbose plotting output (default: false)
 * @return RangeData with calculated ranges r1Tx, r2Tx, r1Rx, r2Rx
 */
RangeData calculaSlantRangeFermat(const ProcessedDEM& demData,
                                 const TrajectoryMatrix& trajectories,
                                 const Environment& env,
                                 ReflectionPoints& strReflexao,
                                 RefractionPoints& strRefraccoes,
                                 bool bPlotVerbose = false);

/**
 * @brief Create CHIRP signal
 * 
 * Generates linear frequency modulated chirp signal
 * Equivalent to MATLAB: chirp = createCHRIP(strRadarTx);
 */
std::vector<std::complex<double>> createCHIRP(const RadarTx& radarTx);

/**
 * @brief Process raw radar data
 * 
 * Complete raw data processing pipeline including:
 * - Time and phase calculations
 * - Raw data matrix generation
 * - Pulse compression
 * - FFT processing
 */
RawDataResult processRawData(const TrajectoryMatrix& trajectories,
                            const RangeData& ranges,
                            const RadarTx& radarTx,
                            const RadarRx& radarRx,
                            const Target& target,
                            const System& system,
                            const Environment& env);

// Helper functions for Fermat calculations
std::vector<double> calculateRefractionPoint(double x1, double y1, double z1,
                                            double x2, double y2, double z2,
                                            const ProcessedDEM& demData,
                                            double n1, double n2,
                                            bool bPlotVerbose = false);

std::vector<double> calculateReflectionPoint(double txX, double txY, double txZ,
                                           double rxX, double rxY, double rxZ,
                                           const ProcessedDEM& demData);

// Additional helper functions for DEM interpolation
double interpolateDEMElevation(double x, double y, const ProcessedDEM& demData);

// Helper function for FFT-based correlation
std::vector<std::complex<double>> correlateWithReference(
    const std::vector<std::complex<double>>& signal,
    const std::vector<std::complex<double>>& reference);

#endif // RAW_DATA_PROCESSING_H