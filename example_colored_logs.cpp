// Ejemplo de uso en main.cpp o cualquier archivo:

#include "color_logger.h"

int main() {
    LOG_INFO("Starting radar trajectory calculation...");
    LOG_PROGRESS("Loading DEM file...");
    LOG_SUCCESS("DEM loaded successfully!");
    LOG_TIMER("Processing time: 125 ms");
    LOG_WARNING("Using debug mode with only 10 points");
    LOG_ERROR("Failed to load configuration file");
    
    return 0;
}