#ifndef LOAD_DEM_H
#define LOAD_DEM_H

#include "radar_params.h"
#include <vector>
#include <string>

/**
 * @brief Load Digital Elevation Model from GeoTIFF file
 * 
 * This function reads a GeoTIFF file and converts the geographic coordinates
 * to UTM coordinates, centering them at the origin.
 * 
 * @param filename Path to the GeoTIFF (.tif) file
 * @return ProcessedDEM structure containing DEM data
 */
ProcessedDEM loadDEM(const std::string& filename);

#endif // LOAD_DEM_H