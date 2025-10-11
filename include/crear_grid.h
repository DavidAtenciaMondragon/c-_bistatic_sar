#ifndef CREAR_GRID_H
#define CREAR_GRID_H

#include "radar_params.h"
#include <vector>

// Structure to hold 3D meshgrid data
struct Grid3D {
    std::vector<std::vector<std::vector<double>>> X; // X coordinates meshgrid
    std::vector<std::vector<std::vector<double>>> Y; // Y coordinates meshgrid  
    std::vector<std::vector<std::vector<double>>> Z; // Z coordinates meshgrid
    GridToProc gridToProc;                           // Axis information
    int nx, ny, nz;                                  // Grid dimensions
};

/**
 * @brief Create 3D processing grid around target position
 * 
 * This function creates a 3D grid centered around the target position
 * equivalent to the MATLAB code:
 * 
 * x0 = strTarget.pos(1);
 * y0 = strTarget.pos(2);  dxy = 0.04;  lxy = 0.8;
 * z0 = strTarget.pos(3);   dz = 0.08;   lz = 0.8;
 * 
 * xAxis = x0-lxy:dxy:x0+lxy;
 * yAxis = y0-lxy:dxy:y0+lxy;  
 * zAxis = z0-lz:dz:z0+lz;
 * 
 * [X,Y,Z] = meshgrid(xAxis,yAxis,zAxis);
 * 
 * @param target Target parameters containing position
 * @param dxy Grid spacing in X and Y directions [m]
 * @param lxy Half-size of grid in X and Y directions [m]
 * @param dz Grid spacing in Z direction [m]
 * @param lz Half-size of grid in Z direction [m]
 * 
 * @return Grid3D structure containing meshgrid and axis data
 */
Grid3D crearGrid(const Target& target,
                 double dxy = 0.04,  // Grid spacing XY [m]
                 double lxy = 0.8,   // Half-size XY [m]
                 double dz = 0.08,   // Grid spacing Z [m]
                 double lz = 0.8);   // Half-size Z [m]

#endif // CREAR_GRID_H