#include "crear_grid.h"
#include "color_logger.h"  // Agregar ColorLogger
#include <iostream>
#include <iomanip>

Grid3D crearGrid(const Target& target, double dxy, double lxy, double dz, double lz) {
    
    std::cout << "Creating 3D processing grid..." << std::endl;
    
    Grid3D grid;
    
    // Extract target position
    double x0 = target.pos[0];
    double y0 = target.pos[1]; 
    double z0 = target.pos[2];
    
    std::cout << "Target position: (" << x0 << ", " << y0 << ", " << z0 << ")" << std::endl;
    std::cout << "Grid parameters: dxy=" << dxy << ", lxy=" << lxy << ", dz=" << dz << ", lz=" << lz << std::endl;
    
    // Create axes vectors (equivalent to MATLAB: x0-lxy:dxy:x0+lxy)
    
    // X axis: x0-lxy:dxy:x0+lxy
    for (double x = x0 - lxy; x <= x0 + lxy + dxy/2; x += dxy) {
        grid.gridToProc.xAxis.push_back(x);
    }
    
    // Y axis: y0-lxy:dxy:y0+lxy  
    for (double y = y0 - lxy; y <= y0 + lxy + dxy/2; y += dxy) {
        grid.gridToProc.yAxis.push_back(y);
    }
    
    // Z axis: z0-lz:dz:z0+lz
    for (double z = z0 - lz; z <= z0 + lz + dz/2; z += dz) {
        grid.gridToProc.zAxis.push_back(z);
    }
    
    // Store dimensions
    grid.nx = grid.gridToProc.xAxis.size();
    grid.ny = grid.gridToProc.yAxis.size(); 
    grid.nz = grid.gridToProc.zAxis.size();
    
    std::cout << "Grid dimensions: " << grid.nx << " x " << grid.ny << " x " << grid.nz 
              << " = " << (grid.nx * grid.ny * grid.nz) << " points" << std::endl;
    
    // Create 3D meshgrid (equivalent to MATLAB: [X,Y,Z] = meshgrid(xAxis,yAxis,zAxis))
    // Note: MATLAB meshgrid has specific dimension ordering
    
    // Initialize 3D vectors
    grid.X.resize(grid.ny, std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.nz)));
    grid.Y.resize(grid.ny, std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.nz)));
    grid.Z.resize(grid.ny, std::vector<std::vector<double>>(grid.nx, std::vector<double>(grid.nz)));
    
    // Fill meshgrid (MATLAB meshgrid convention: [Y, X, Z] indexing)
    for (int j = 0; j < grid.ny; j++) {        // Y dimension (rows)
        for (int i = 0; i < grid.nx; i++) {    // X dimension (columns)
            for (int k = 0; k < grid.nz; k++) { // Z dimension (pages)
                grid.X[j][i][k] = grid.gridToProc.xAxis[i];  // X varies with column index
                grid.Y[j][i][k] = grid.gridToProc.yAxis[j];  // Y varies with row index
                grid.Z[j][i][k] = grid.gridToProc.zAxis[k];  // Z varies with page index
            }
        }
    }
    
    // Display grid range information
    std::cout << "X range: [" << std::fixed << std::setprecision(3) 
              << grid.gridToProc.xAxis.front() << ", " << grid.gridToProc.xAxis.back() << "]" << std::endl;
    std::cout << "Y range: [" << grid.gridToProc.yAxis.front() << ", " << grid.gridToProc.yAxis.back() << "]" << std::endl;
    std::cout << "Z range: [" << grid.gridToProc.zAxis.front() << ", " << grid.gridToProc.zAxis.back() << "]" << std::endl;
    
    LOG_SUCCESS("+ 3D grid creation completed!");
    
    return grid;
}