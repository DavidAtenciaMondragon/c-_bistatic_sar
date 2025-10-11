#include "load_dem.h"
#include "radar_params.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

ProcessedDEM loadDEM(const std::string& filename) {
    ProcessedDEM dem;
    
    std::cout << "Loading DEM from file: " << filename << std::endl;
    std::cout << "Using simplified DEM loading (no GDAL support)" << std::endl;
    
    // Try to load as simple text format first
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open DEM file: " << filename << std::endl;
        std::cerr << "Creating default DEM data..." << std::endl;
        
        // Create a comprehensive DEM covering all quadrants
        // Parameters: 32x32 grid from -500 to +500 in both X and Y
        const double lim = 500.0;  // Spatial limit in meters
        const int M = 32;          // Grid size (MxM matrix)
        
        int ncols = M;
        int nrows = M;
        double xmin = -lim;
        double xmax = +lim;
        double ymin = -lim;
        double ymax = +lim;
        double cellsize = (xmax - xmin) / (ncols - 1);  // Grid spacing
        
        dem.Z_DEM.resize(nrows, std::vector<double>(ncols));
        
        // Create X_vec and Y_vec covering all quadrants
        dem.X_vec.resize(ncols);
        dem.Y_vec.resize(nrows);
        
        // X coordinates: from -500 to +500
        for (int i = 0; i < ncols; i++) {
            dem.X_vec[i] = xmin + i * cellsize;
        }
        
        // Y coordinates: from -500 to +500  
        for (int i = 0; i < nrows; i++) {
            dem.Y_vec[i] = ymin + i * cellsize;
        }
        
        // Create terrain with elevation variation covering all quadrants
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                double x = dem.X_vec[j];
                double y = dem.Y_vec[i];
                
                // Create interesting terrain with multiple features
                // Base elevation with gentle slope
                double base_elevation = 0.0;
                
                // Add sinusoidal variation
                double wave1 = 0.0 * sin(x * 0.01) * cos(y * 0.01);
                double wave2 = 0.0 * sin(x * 0.02) * sin(y * 0.015);
                
                // Add radial feature (hill/depression based on distance from origin)
                double r = sqrt(x*x + y*y);
                double radial = 0.0 * exp(-r*r / (200.0*200.0));
                
                // Combine all elevation components
                dem.Z_DEM[i][j] = base_elevation + wave1 + wave2 + radial;
            }
        }
        
        std::cout << "Created comprehensive DEM: " << ncols << "x" << nrows << " cells" << std::endl;
        std::cout << "  Coverage: X=[" << xmin << ", " << xmax << "], Y=[" << ymin << ", " << ymax << "]" << std::endl;
        std::cout << "  Cell size: " << cellsize << " meters" << std::endl;
        return dem;
    }
    
    std::string line;
    bool header_complete = false;
    int row = 0;
    int ncols = 0, nrows = 0;
    double xllcorner = 0.0, yllcorner = 0.0, cellsize = 1.0;
    double nodata_value = -9999.0;
    
    // Try to parse as ESRI ASCII Grid format
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        
        if (!header_complete) {
            iss >> token;
            
            if (token == "ncols") {
                iss >> ncols;
            } else if (token == "nrows") {
                iss >> nrows;
            } else if (token == "xllcorner") {
                iss >> xllcorner;
            } else if (token == "yllcorner") {
                iss >> yllcorner;
            } else if (token == "cellsize") {
                iss >> cellsize;
            } else if (token == "NODATA_value" || token == "nodata_value") {
                iss >> nodata_value;
            } else {
                // This line contains data, start reading grid
                header_complete = true;
                dem.Z_DEM.resize(nrows, std::vector<double>(ncols));
                
                // Create X_vec and Y_vec
                dem.X_vec.resize(ncols);
                dem.Y_vec.resize(nrows);
                
                for (int i = 0; i < ncols; i++) {
                    dem.X_vec[i] = xllcorner + i * cellsize;
                }
                
                for (int i = 0; i < nrows; i++) {
                    dem.Y_vec[i] = yllcorner + i * cellsize;
                }
                
                // Parse this first data line
                std::istringstream data_iss(line);
                for (int col = 0; col < ncols; col++) {
                    if (data_iss >> dem.Z_DEM[row][col]) {
                        // Successfully read value
                    } else {
                        dem.Z_DEM[row][col] = nodata_value;
                    }
                }
                row++;
            }
        } else {
            // Reading data rows
            if (row < nrows) {
                for (int col = 0; col < ncols; col++) {
                    if (iss >> dem.Z_DEM[row][col]) {
                        // Successfully read value
                    } else {
                        dem.Z_DEM[row][col] = nodata_value;
                    }
                }
                row++;
            }
        }
    }
    
    file.close();
    
    // Validate that we read the expected amount of data
    if (row != nrows) {
        std::cerr << "Warning: Expected " << nrows << " rows, but read " << row << " rows" << std::endl;
    }
    
    std::cout << "Successfully loaded DEM:" << std::endl;
    std::cout << "  Size: " << ncols << " x " << nrows << " cells" << std::endl;
    std::cout << "  Cell size: " << cellsize << std::endl;
    std::cout << "  Lower-left corner: (" << xllcorner << ", " << yllcorner << ")" << std::endl;
    std::cout << "  NoData value: " << nodata_value << std::endl;
    
    return dem;
}