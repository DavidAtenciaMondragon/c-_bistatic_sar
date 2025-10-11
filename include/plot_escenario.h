#ifndef PLOT_ESCENARIO_H
#define PLOT_ESCENARIO_H

#include "funcao_espiral.h"
#include "load_dem.h"
#include "radar_params.h"
#include <string>

/**
 * @brief Plot scenario with radar trajectories, DEM surface, and target position
 * 
 * This function creates output files with data for visualization of:
 * - Transmitter trajectory (black line in MATLAB)
 * - Receiver trajectory (red line in MATLAB) 
 * - Target position (black marker in MATLAB)
 * - DEM surface (with transparency in MATLAB)
 * 
 * The function generates CSV files that can be loaded into plotting tools
 * like Python matplotlib, MATLAB, or other visualization software.
 * 
 * @param trajTx Transmitter trajectory points
 * @param trajRx Receiver trajectory points  
 * @param demData DEM grid data (X, Y, Z coordinates)
 * @param target Target parameters including position
 * @param outputDir Directory to save output files (default: "plot_data/")
 */
void plotEscenario(const TrajectoryPoints& trajTx, 
                   const TrajectoryPoints& trajRx,
                   const ProcessedDEM& demData,
                   const Target& target,
                   const std::string& outputDir = "plot_data/");

/**
 * @brief Generate Python script to visualize the scenario
 * 
 * Creates a Python script that loads the CSV data and creates
 * a 3D plot equivalent to the MATLAB visualization.
 * 
 * @param outputDir Directory containing the data files
 */
void generatePythonPlotScript(const std::string& outputDir = "plot_data/");

#endif // PLOT_ESCENARIO_H