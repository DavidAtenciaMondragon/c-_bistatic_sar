#include "plot_escenario.h"
#include "color_logger.h"  // Agregar ColorLogger
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

void plotEscenario(const TrajectoryPoints& trajTx, 
                   const TrajectoryPoints& trajRx,
                   const ProcessedDEM& demData,
                   const Target& target,
                   const std::string& outputDir) {
    
    std::cout << "Generating plot data files..." << std::endl;
    
    // Create output directory if it doesn't exist
    mkdir(outputDir.c_str(), 0755);
    
    // 1. Save transmitter trajectory
    std::ofstream txFile(outputDir + "/transmitter_trajectory.csv");
    if (txFile.is_open()) {
        txFile << "X,Y,Z\n";
        for (size_t i = 0; i < trajTx.Px.size(); i++) {
            txFile << std::fixed << std::setprecision(6)
                   << trajTx.Px[i] << "," 
                   << trajTx.Py[i] << "," 
                   << trajTx.Pz[i] << "\n";
        }
        txFile.close();
        LOG_SUCCESS("+ Transmitter trajectory saved (" + std::to_string(trajTx.Px.size()) + " points)");
    }
    
    // 2. Save receiver trajectory  
    std::ofstream rxFile(outputDir + "/receiver_trajectory.csv");
    if (rxFile.is_open()) {
        rxFile << "X,Y,Z\n";
        for (size_t i = 0; i < trajRx.Px.size(); i++) {
            rxFile << std::fixed << std::setprecision(6)
                   << trajRx.Px[i] << "," 
                   << trajRx.Py[i] << "," 
                   << trajRx.Pz[i] << "\n";
        }
        rxFile.close();
        LOG_SUCCESS("+ Receiver trajectory saved (" + std::to_string(trajRx.Px.size()) + " points)");
    }
    
    // 3. Save target position
    std::ofstream targetFile(outputDir + "/target_position.csv");
    if (targetFile.is_open()) {
        targetFile << "X,Y,Z\n";
        targetFile << std::fixed << std::setprecision(6)
                   << target.pos[0] << "," 
                   << target.pos[1] << "," 
                   << target.pos[2] << "\n";
        targetFile.close();
        LOG_SUCCESS("+ Target position saved");
    }
    
    // 4. Save DEM surface data
    std::ofstream demFile(outputDir + "/dem_surface.csv");
    if (demFile.is_open() && !demData.X_vec.empty() && !demData.Y_vec.empty()) {
        demFile << "X,Y,Z\n";
        for (size_t row = 0; row < demData.Y_vec.size(); row++) {
            for (size_t col = 0; col < demData.X_vec.size(); col++) {
                if (row < demData.Z_DEM.size() && col < demData.Z_DEM[row].size()) {
                    demFile << std::fixed << std::setprecision(6)
                            << demData.X_vec[col] << "," 
                            << demData.Y_vec[row] << "," 
                            << demData.Z_DEM[row][col] << "\n";
                }
            }
        }
        demFile.close();
        LOG_SUCCESS("+ DEM surface saved (" + std::to_string(demData.Y_vec.size()) + "x" + std::to_string(demData.X_vec.size()) + " points)");
    }
    
    // 5. Save DEM grid structure for surface plotting
    std::ofstream demGridFile(outputDir + "/dem_grid_info.txt");
    if (demGridFile.is_open() && !demData.X_vec.empty() && !demData.Y_vec.empty()) {
        demGridFile << "nRows=" << demData.Y_vec.size() << "\n";
        demGridFile << "nCols=" << demData.X_vec.size() << "\n";
        demGridFile.close();
    }
    
    // Generate Python plotting script
    generatePythonPlotScript(outputDir);
    
    LOG_SUCCESS("+ Plot data generation completed!");
    std::cout << "Files saved in: " << outputDir << std::endl;
    std::cout << "Run 'python plot_scenario.py' to visualize the data" << std::endl;
}

void generatePythonPlotScript(const std::string& outputDir) {
    std::ofstream scriptFile(outputDir + "/plot_scenario.py");
    if (scriptFile.is_open()) {
        scriptFile << R"(#!/usr/bin/env python3
"""
Radar Scenario Visualization Script
Generated automatically from C++ plotEscenario function

This script replicates the MATLAB plot:
figure
hold on 
plot3(PxT, PyT, PzT,'k')
plot3(PxR, PyR, PzR,'r')
plot3(strTarget.pos(1),strTarget.pos(2),strTarget.pos(3),'o','MarkerFaceColor','k');
surf(X_DEM,Y_DEM,Z_DEM,'EdgeColor','none','FaceAlpha',0.1)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid minor 
hold off
legend("Tx","Rx","Target")
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read data files
try:
    tx_data = pd.read_csv('transmitter_trajectory.csv')
    rx_data = pd.read_csv('receiver_trajectory.csv')
    target_data = pd.read_csv('target_position.csv')
    dem_data = pd.read_csv('dem_surface.csv')
    
    # Read grid info
    with open('dem_grid_info.txt', 'r') as f:
        lines = f.readlines()
        nRows = int(lines[0].split('=')[1])
        nCols = int(lines[1].split('=')[1])
    
    print(f"Data loaded successfully!")
    print(f"Transmitter points: {len(tx_data)}")
    print(f"Receiver points: {len(rx_data)}")
    print(f"DEM grid: {nRows}x{nCols}")
    
except Exception as e:
    print(f"Error loading data: {e}")
    exit(1)

# Create figure
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# Plot transmitter trajectory (black line)
ax.plot(tx_data['X'], tx_data['Y'], tx_data['Z'], 'k-', linewidth=1.5, label='Tx')

# Plot receiver trajectory (red line)  
ax.plot(rx_data['X'], rx_data['Y'], rx_data['Z'], 'r-', linewidth=1.5, label='Rx')

# Plot target position (black marker)
ax.scatter(target_data['X'], target_data['Y'], target_data['Z'], 
          c='black', s=100, marker='o', label='Target')

# Plot DEM surface
X_dem = dem_data['X'].values.reshape(nRows, nCols)
Y_dem = dem_data['Y'].values.reshape(nRows, nCols)  
Z_dem = dem_data['Z'].values.reshape(nRows, nCols)

ax.plot_surface(X_dem, Y_dem, Z_dem, alpha=0.1, cmap='terrain', 
               linewidth=0, antialiased=True)

# Set labels and formatting
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.grid(True, alpha=0.3)
ax.legend()

# Set title
ax.set_title('Radar Scenario: Transmitter/Receiver Trajectories with DEM')

# Show plot
plt.tight_layout()
plt.show()

print("Plot displayed successfully!")
)";
        scriptFile.close();
        LOG_SUCCESS("+ Python plotting script generated");
    }
}