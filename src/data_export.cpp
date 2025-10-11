#include "data_export.h"
#include "color_logger.h"  // Agregar ColorLogger
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <cstdint>    // For uint32_t, uint64_t
#include <cstring>    // For strncpy
#include <algorithm>  // For std::min

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

void exportProcessedData(const std::string& outputDir,
                        const RawDataResult& rawResult,
                        const TrajectoryMatrix& trajectories,
                        const ProcessedDEM& demData,
                        const Grid3D& grid,
                        const Environment& env) {
    
    std::cout << "Exporting processed data..." << std::endl;
    
    // Create output directory
    std::string procDir = outputDir + "/proc";
    mkdir(outputDir.c_str(), 0755);
    mkdir(procDir.c_str(), 0755);
    
    // Export raw data results as complex data for processing
    // Raw data (uncompressed): Export complex values
    exportComplexMatrix(procDir + "/rawData_complex.bin", rawResult.rawData);
    
    // Root data (compressed): Export complex values
    exportComplexMatrix(procDir + "/rootData_complex.bin", rawResult.rootData);
    
    // Also export auxiliary data as complex values
    exportComplexMatrix(procDir + "/auxData_complex.bin", rawResult.auxData);
    
    // Export trajectory data
    exportRealMatrix(procDir + "/Tx_positions.csv", trajectories.Tx);
    exportRealMatrix(procDir + "/Rx_positions.csv", trajectories.Rx);
    exportVector(procDir + "/target_position.csv", trajectories.P);
    
    // Export DEM data
    if (!demData.X_vec.empty()) {
        exportVector(procDir + "/DEM_X_vec.csv", demData.X_vec);
        exportVector(procDir + "/DEM_Y_vec.csv", demData.Y_vec);
        exportRealMatrix(procDir + "/DEM_Z_matrix.csv", demData.Z_DEM);
    }
    
    // Export grid data
    exportVector(procDir + "/grid_xAxis.csv", grid.gridToProc.xAxis);
    exportVector(procDir + "/grid_yAxis.csv", grid.gridToProc.yAxis);
    exportVector(procDir + "/grid_zAxis.csv", grid.gridToProc.zAxis);
    
    // Export environment parameters
    std::ofstream envFile(procDir + "/environment.txt");
    if (envFile.is_open()) {
        envFile << "n1=" << env.n1 << std::endl;
        envFile << "n2=" << env.n2 << std::endl;
        envFile.close();
    }
    
    // Generate Python visualization script - DISABLED (using independent scripts)
    // generateVisualizationScript(procDir);
    
    LOG_SUCCESS("+ Processed data exported to: " + procDir);
}

void exportComplexMatrix(const std::string& filename,
                        const std::vector<std::vector<std::complex<double>>>& data) {
    if (data.empty() || data[0].empty()) {
        std::cerr << "Error: Empty data matrix" << std::endl;
        return;
    }
    
    // Open binary file for writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    try {
        // Matrix dimensions
        uint32_t numDimensions = 2;  // 2D matrix
        uint64_t rows = static_cast<uint64_t>(data.size());
        uint64_t cols = static_cast<uint64_t>(data[0].size());
        
        // Header information
        uint32_t isComplex = 1;  // Complex data
        std::string dataTypeStr = "double";
        uint32_t dataTypeLen = static_cast<uint32_t>(dataTypeStr.length());
        
        // Write header
        // 1. Number of dimensions (uint32)
        file.write(reinterpret_cast<const char*>(&numDimensions), sizeof(uint32_t));
        
        // 2. Complexity flag (uint32): 0=Real, 1=Complex
        file.write(reinterpret_cast<const char*>(&isComplex), sizeof(uint32_t));
        
        // 3. Data type string length (uint32)
        file.write(reinterpret_cast<const char*>(&dataTypeLen), sizeof(uint32_t));
        
        // 4. Data type string (char[16] fixed, padded with nulls)
        char fixedDataType[16] = {0};  // Initialize with nulls
        strncpy(fixedDataType, dataTypeStr.c_str(), std::min(dataTypeStr.length(), size_t(15)));
        file.write(fixedDataType, 16);
        
        // 5. Write dimensions (uint64[])
        file.write(reinterpret_cast<const char*>(&rows), sizeof(uint64_t));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(uint64_t));
        
        // 6. Write data (interleaved Real-Imaginary for complex)
        for (size_t i = 0; i < data.size(); i++) {
            for (size_t j = 0; j < data[i].size(); j++) {
                double realPart = data[i][j].real();
                double imagPart = data[i][j].imag();
                // Write real part
                file.write(reinterpret_cast<const char*>(&realPart), sizeof(double));
                // Write imaginary part
                file.write(reinterpret_cast<const char*>(&imagPart), sizeof(double));
            }
        }
        
        LOG_SUCCESS("✅ Complex matrix saved to binary: " + filename);
        LOG_INFO("   Dimensions: " + std::to_string(rows) + "x" + std::to_string(cols));
        LOG_INFO("   Type: double (Complex, R-I interleaved)");
        
    } catch (const std::exception& e) {
        std::cerr << "Error writing binary file: " << e.what() << std::endl;
        file.close();
        return;
    }
    
    file.close();
}

void exportRealMatrix(const std::string& filename,
                     const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << std::fixed << std::setprecision(6);
    
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < data[i].size(); j++) {
            if (j > 0) file << ",";
            file << data[i][j];
        }
        file << std::endl;
    }
    
    file.close();
}

void exportVector(const std::string& filename,
                 const std::vector<double>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << std::fixed << std::setprecision(6);
    
    for (size_t i = 0; i < data.size(); i++) {
        if (i > 0) file << ",";
        file << data[i];
    }
    file << std::endl;
    
    file.close();
}

void generateProcessingReport(const std::string& outputDir,
                             const RadarTx& radarTx,
                             const RadarRx& radarRx,
                             const Target& target,
                             const System& system,
                             const Environment& env,
                             const TrajectoryMatrix& trajectories,
                             const Grid3D& grid) {
    
    std::string reportFile = outputDir + "/processing_report.txt";
    std::ofstream report(reportFile);
    
    if (!report.is_open()) {
        std::cerr << "Error: Could not create report file" << std::endl;
        return;
    }
    
    report << "=== RADAR PROCESSING REPORT ===" << std::endl;
    report << std::fixed << std::setprecision(3);
    
    // System parameters
    report << "\nSYSTEM PARAMETERS:" << std::endl;
    report << "Speed of light: " << system.VelocidadeLuz/1e8 << " x 10^8 m/s" << std::endl;
    report << "Maximum index: " << system.IndiceMaximo << std::endl;
    report << "Number of pulses: " << system.NumeroPulsos << std::endl;
    
    // Radar parameters
    report << "\nRADAR PARAMETERS:" << std::endl;
    report << "Carrier frequency: " << radarTx.FreqPortadora/1e9 << " GHz" << std::endl;
    report << "Wavelength: " << radarTx.lamb*1000 << " mm" << std::endl;
    report << "PRF: " << radarTx.PRF << " Hz" << std::endl;
    report << "Sampling frequency: " << radarTx.fs/1e6 << " MHz" << std::endl;
    report << "Pulse duration: " << radarTx.DuracaoPulso*1e6 << " μs" << std::endl;
    
    // Target parameters
    report << "\nTARGET PARAMETERS:" << std::endl;
    report << "Position: (" << target.pos[0] << ", " << target.pos[1] << ", " << target.pos[2] << ") m" << std::endl;
    report << "RCS: " << target.rcs << " m²" << std::endl;
    
    // Environment
    report << "\nENVIRONMENT:" << std::endl;
    report << "Refractive index 1 (air): " << env.n1 << std::endl;
    report << "Refractive index 2 (ground): " << env.n2 << std::endl;
    
    // Trajectory statistics
    report << "\nTRAJECTORY STATISTICS:" << std::endl;
    report << "Number of positions: " << trajectories.Tx.size() << std::endl;
    
    if (!trajectories.Tx.empty()) {
        double minZ = trajectories.Tx[0][2], maxZ = trajectories.Tx[0][2];
        for (const auto& pos : trajectories.Tx) {
            minZ = std::min(minZ, pos[2]);
            maxZ = std::max(maxZ, pos[2]);
        }
        report << "Altitude range: [" << minZ << ", " << maxZ << "] m" << std::endl;
    }
    
    // Grid statistics
    report << "\nPROCESSING GRID:" << std::endl;
    report << "Grid dimensions: " << grid.nx << " x " << grid.ny << " x " << grid.nz << std::endl;
    report << "Total grid points: " << (grid.nx * grid.ny * grid.nz) << std::endl;
    
    if (!grid.gridToProc.xAxis.empty()) {
        report << "X range: [" << grid.gridToProc.xAxis.front() << ", " << grid.gridToProc.xAxis.back() << "] m" << std::endl;
        report << "Y range: [" << grid.gridToProc.yAxis.front() << ", " << grid.gridToProc.yAxis.back() << "] m" << std::endl;
        report << "Z range: [" << grid.gridToProc.zAxis.front() << ", " << grid.gridToProc.zAxis.back() << "] m" << std::endl;
    }
    
    report << "\n=== END OF REPORT ===" << std::endl;
    report.close();
    
    LOG_SUCCESS("+ Processing report generated: " + reportFile);
}
