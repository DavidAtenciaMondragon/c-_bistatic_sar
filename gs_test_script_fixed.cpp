#include <iostream>
#include <cmath>
#include <vector>
#include "radar_params.h"
#include "radar_init.h"
#include "trajectory_generator.h"
#include "load_dem.h"
#include "plot_escenario.h"
#include "crear_grid.h"
#include "raw_data_processing.h"
#include "data_export.h"
#include "color_logger.h"  // Agregar ColorLogger

int main() {
    // Clear screen equivalent
    LOG_INFO("> Starting radar trajectory calculation...");
    
    // Initialize parameters 
    RadarTx strRadarTx = initializeRadarTx();
    RadarRx strRadarRx = initializeRadarRx();
    
    // Initialize target with default position (origin)
    Target strTarget = initializeTarget();
    
    // For gs_test_script, set a specific target position
    // TESTING: Try with intermediate position to get closer to rngBin=85
    
    
    System strSystem = initializeSystem();
    
    // Calculate wavelength
    strRadarTx.lamb = strSystem.VelocidadeLuz / strRadarTx.FreqPortadora;
    
    std::cout << "Calculated wavelength: " << strRadarTx.lamb << " m" << std::endl;
    
    // Create trajectory section
    std::cout << "Creating trajectories..." << std::endl;
    
    // Generate transmitter and receiver trajectories
    TrajectoryPoints trajTx = generateTransmitterTrajectory(strRadarTx);
    TrajectoryPoints trajRx = generateReceiverTrajectory(strRadarRx);

    // For debug, use just the first 10 points of trajectories
    // trajTx.Px.resize(10);
    // trajTx.Py.resize(10);
    // trajTx.Pz.resize(10);
    // trajRx.Px.resize(10);
    // trajRx.Py.resize(10);
    // trajRx.Pz.resize(10);
    
    // Display some results
    std::cout << "Transmitter trajectory points: " << trajTx.Px.size() << std::endl;
    std::cout << "Receiver trajectory points: " << trajRx.Px.size() << std::endl;
    
    // Load .tif file section
    std::cout << "\nLoading DEM file..." << std::endl;
    std::string name_DEM = "assets/DEM_1x1km_Res30m_Lat-3_6160_Lon-80_4552.tif";
    
    ProcessedDEM demData = loadDEM(name_DEM);
    
    if (demData.Z_DEM.size() > 0 && demData.X_vec.size() > 0) {
        std::cout << "DEM loaded successfully!" << std::endl;
        std::cout << "DEM dimensions: " << demData.Y_vec.size() << " x " << demData.X_vec.size() << std::endl;
        
        // Display some DEM statistics
        double minZ = 1e9, maxZ = -1e9;
        for (size_t row = 0; row < demData.Z_DEM.size(); row++) {
            for (size_t col = 0; col < demData.Z_DEM[row].size(); col++) {
                minZ = std::min(minZ, demData.Z_DEM[row][col]);
                maxZ = std::max(maxZ, demData.Z_DEM[row][col]);
            }
        }
        std::cout << "Elevation range: " << minZ << " to " << maxZ << " meters" << std::endl;
    } else {
        std::cout << "Warning: DEM file could not be loaded. Check file path and GDAL installation." << std::endl;
    }
    
    // Plot scenario section
    std::cout << "\nGenerating scenario plot..." << std::endl;
    plotEscenario(trajTx, trajRx, demData, strTarget);
    
    // Create grid section
    std::cout << "\nCreating processing grid..." << std::endl;
    
    // Grid parameters (from MATLAB code)
    double dxy = 0.04;  // Grid spacing in X and Y [m]
    double lxy = 0.8;   // Half-size in X and Y [m] 
    double dz = 0.08;   // Grid spacing in Z [m]
    double lz = 0.8;    // Half-size in Z [m]
    
    Grid3D grid = crearGrid(strTarget, dxy, lxy, dz, lz);
    
    // Display some grid statistics
    std::cout << "Grid center point: (" 
              << grid.X[grid.ny/2][grid.nx/2][grid.nz/2] << ", "
              << grid.Y[grid.ny/2][grid.nx/2][grid.nz/2] << ", "
              << grid.Z[grid.ny/2][grid.nx/2][grid.nz/2] << ")" << std::endl;
    
    // Raw data processing section
    std::cout << "\n=== RAW DATA PROCESSING ===" << std::endl;
    
    // Initialize environment parameters
    Environment strEnvironment = initializeEnvironment();
    
    // Prepare trajectory data for processing
    TrajectoryMatrix trajectories = prepareTrajectoryData(trajTx, trajRx, strTarget);
    
    // Prepare DEM data for processing  
    ProcessedDEM processedDEM = prepareDEMData(demData);
    
    // Calculate ranges using Fermat's principle for refraction
    ReflectionPoints strReflexao;
    RefractionPoints strRefraccoes; 
    RangeData ranges = calculaSlantRangeFermat(processedDEM, trajectories, strEnvironment, 
                                              strReflexao, strRefraccoes, false);
    
    // Process raw radar data
    RawDataResult rawResult = processRawData(trajectories, ranges, strRadarTx, strRadarRx, 
                                           strTarget, strSystem, strEnvironment);
    
    // Export all processed data
    std::cout << "\n=== DATA EXPORT ===" << std::endl;
    exportProcessedData("output", rawResult, trajectories, processedDEM, grid, strEnvironment);
    
    // Generate processing report
    generateProcessingReport("output", strRadarTx, strRadarRx, strTarget, strSystem, 
                           strEnvironment, trajectories, grid);
    
    std::cout << "\n[SUCCESS] All radar processing completed successfully!" << std::endl;
    LOG_SUCCESS("+ All radar processing completed successfully!");
    std::cout << "Check 'output/proc/' folder for processed data files." << std::endl;
    
    return 0;
}
