#include "load_dem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

// Versión simplificada sin GDAL - solo para demostración
DEMData loadDEM(const std::string& filename) {
    DEMData demData;
    
    std::cout << "Note: This is a simplified version of loadDEM without GDAL." << std::endl;
    std::cout << "To use real GeoTIFF files, install GDAL library." << std::endl;
    std::cout << "Attempted to load: " << filename << std::endl;
    
    // Crear datos de ejemplo para demostración
    int nRows = 100;
    int nCols = 100;
    
    demData.nRows = nRows;
    demData.nCols = nCols;
    
    // Inicializar grillas
    demData.Xg.resize(nRows, std::vector<double>(nCols));
    demData.Yg.resize(nRows, std::vector<double>(nCols));
    demData.Zg.resize(nRows, std::vector<double>(nCols));
    
    // Generar datos de ejemplo (terreno sintético)
    double range = 1000.0; // rango de 1km
    double dx = range / (nCols - 1);
    double dy = range / (nRows - 1);
    
    for (int row = 0; row < nRows; row++) {
        for (int col = 0; col < nCols; col++) {
            // Coordenadas centradas en origen
            demData.Xg[row][col] = -range/2.0 + col * dx;
            demData.Yg[row][col] = -range/2.0 + row * dy;
            
            // Elevación sintética (colina suave)
            double x = demData.Xg[row][col];
            double y = demData.Yg[row][col];
            double r = sqrt(x*x + y*y);
            demData.Zg[row][col] = 0.0 * exp(-r*r / (2.0 * 300.0 * 300.0));
        }
    }
    
    std::cout << "Generated synthetic DEM: " << nRows << "x" << nCols << " points" << std::endl;
    
    return demData;
}