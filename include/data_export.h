#ifndef DATA_EXPORT_H
#define DATA_EXPORT_H

#include "raw_data_processing.h"
#include "crear_grid.h"
#include <string>

/**
 * @brief Export all processed data to files
 * 
 * Creates organized output files equivalent to MATLAB save commands:
 * save(dir, 'rootData', 'strDEM', 'strEnvironment', 'Tx', 'Rx', 'strGridToProc');
 * 
 * @param outputDir Directory to save files
 * @param rawResult Raw data processing results
 * @param trajectories Trajectory matrices
 * @param demData Processed DEM data
 * @param grid 3D grid data
 * @param env Environment parameters
 */
void exportProcessedData(const std::string& outputDir,
                        const RawDataResult& rawResult,
                        const TrajectoryMatrix& trajectories,
                        const ProcessedDEM& demData,
                        const Grid3D& grid,
                        const Environment& env);

/**
 * @brief Export complex matrix to CSV file
 * 
 * Saves complex data with separate real and imaginary parts
 */
void exportComplexMatrix(const std::string& filename,
                        const std::vector<std::vector<std::complex<double>>>& data);

/**
 * @brief Export real matrix to CSV file
 */
void exportRealMatrix(const std::string& filename,
                     const std::vector<std::vector<double>>& data);

/**
 * @brief Export vector to CSV file
 */
void exportVector(const std::string& filename,
                 const std::vector<double>& data);

/**
 * @brief Export 3D complex matrix (SAR output) to binary file
 * 
 * Saves the 3D complex matrix from back projection processing
 * Format: Binary file with header containing dimensions and complex data
 * 
 * @param filename Output binary file path
 * @param output3D 3D matrix with dimensions [nx][ny][nz] containing complex values
 * @param nx X dimension size
 * @param ny Y dimension size  
 * @param nz Z dimension size
 * @param gridInfo Optional grid coordinate information for metadata
 */
void export3DComplexMatrix(const std::string& filename,
                          const std::vector<std::vector<std::vector<std::complex<double>>>>& output3D,
                          size_t nx, size_t ny, size_t nz,
                          const std::vector<double>& grid_x = {},
                          const std::vector<double>& grid_y = {},
                          const std::vector<double>& grid_z = {});

/**
 * @brief Generate processing summary report
 */
void generateProcessingReport(const std::string& outputDir,
                             const RadarTx& radarTx,
                             const RadarRx& radarRx,
                             const Target& target,
                             const System& system,
                             const Environment& env,
                             const TrajectoryMatrix& trajectories,
                             const Grid3D& grid);

/**
 * @brief Generate Python visualization script for radar data
 */

#endif // DATA_EXPORT_H"