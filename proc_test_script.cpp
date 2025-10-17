#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <cstdio>  // Para printf
#include <cmath>   // Para M_PI y funciones matemáticas
#include <algorithm>  // Para std::min_element, std::max_element
#define _USE_MATH_DEFINES  // Para M_PI en Windows
#include "radar_params.h"
#include "radar_init.h"
#include "trajectory_generator.h"
#include "binary_reader.h"
#include "interpolation.h" 
#include "csv_reader.h"
#include "matrix_reader.h"
#include "raw_data_processing.h"
#include "color_logger.h"
#include "data_export.h" // Para exportar matriz 3D output
#include "project_paths.h" // Centralized paths
#include <omp.h>  // Para paralelización OpenMP

// Definir M_PI si no está disponible
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Color logger para output mejorado
// Comentado para usar el ColorLogger con colores apropiados
// #define LOG_INFO(msg) std::cout << "[INFO] " << msg << std::endl
// #define LOG_ERROR(msg) std::cerr << "[ERROR] " << msg << std::endl  
// #define LOG_PROGRESS(msg) std::cout << "[PROGRESS] " << msg << std::endl
#define LOG_PERF(msg) std::cout << "[PERFORMANCE] " << msg << std::endl

// ============================================================================
// CONFIGURACION DE VERBOSIDAD
// ============================================================================

struct VerbosityFlags {
    bool show_detailed_logs = false;      // Logs detallados de cálculos
    bool show_numerical_results = false;  // Resultados numéricos específicos
    bool show_debug_info = false;         // Debug de interpolación y acceso 2D
    bool show_refraction_points = false;  // Puntos de refracción específicos
    bool show_coordinate_mapping = false; // Mapeo de coordenadas 3D
    
    // Configuraciones predefinidas
    static VerbosityFlags minimal() {
        VerbosityFlags flags;
        // Solo logs esenciales de progreso
        return flags;
    }
    
    static VerbosityFlags standard() {
        VerbosityFlags flags;
        flags.show_numerical_results = true;
        return flags;
    }
    
    static VerbosityFlags debug() {
        VerbosityFlags flags;
        flags.show_detailed_logs = true;
        flags.show_numerical_results = true;
        flags.show_debug_info = true;
        flags.show_refraction_points = true;
        flags.show_coordinate_mapping = true;
        return flags;
    }
};

// ============================================================================
// VECTORIZED MATHEMATICAL FUNCTIONS (MATLAB-style optimization)
// ============================================================================


/**
 * @brief Vectorized distance calculation for multiple points
 * Equivalent to MATLAB: sqrt(sum((A - B).^2, 2)) for multiple points
 */
inline std::vector<double> vectorizedDistance3D(
    const std::vector<std::vector<double>>& points1,
    const std::vector<std::vector<double>>& points2) {
    
    std::vector<double> distances;
    distances.reserve(std::min(points1.size(), points2.size()));
    
    size_t n = std::min(points1.size(), points2.size());
    
    // SIMD-friendly loop: process multiple distances simultaneously
    for (size_t i = 0; i < n; i++) {
        const auto& p1 = points1[i];
        const auto& p2 = points2[i];
        
        // Vectorized difference calculation
        double dx = p1[0] - p2[0];
        double dy = p1[1] - p2[1];
        double dz = p1[2] - p2[2];
        
        // Fast square root using optimized calculation
        double dist_squared = dx*dx + dy*dy + dz*dz;
        distances.push_back(std::sqrt(dist_squared));
    }
    
    return distances;
}

/**
 * @brief Vectorized optical path calculation for multiple targets
 * Equivalent to MATLAB: n1*d1 + n2*d2 for arrays
 */
inline std::vector<double> vectorizedOpticalPath(
    const std::vector<double>& distances1,
    const std::vector<double>& distances2,
    double n1, double n2) {
    
    std::vector<double> optical_paths;
    optical_paths.reserve(distances1.size());
    
    // SIMD-friendly vectorized calculation
    for (size_t i = 0; i < distances1.size() && i < distances2.size(); i++) {
        optical_paths.push_back(n1 * distances1[i] + n2 * distances2[i]);
    }
    
    return optical_paths;
}

/**
 * @brief Fast grid coordinate generation (MATLAB meshgrid equivalent)
 * More cache-friendly than nested loops
 */
inline void generateFlatGridCoordinates(
    const std::vector<double>& grid_x,
    const std::vector<double>& grid_y, 
    const std::vector<double>& grid_z,
    std::vector<std::vector<double>>& flat_coordinates) {
    
    size_t total_points = grid_x.size() * grid_y.size() * grid_z.size();
    flat_coordinates.clear();
    flat_coordinates.reserve(total_points);
    
    // Generate coordinates in memory-efficient order
    // Changed order: Z -> X -> Y (so Y changes first, then X, then Z)
    for (const auto& z : grid_z) {
        for (const auto& x : grid_x) {
            for (const auto& y : grid_y) {
                flat_coordinates.push_back({x, y, z});
            }
        }
    }
}

/**
 * @brief Batch processing helper for memory efficiency
 * Processes multiple grid points with shared trajectory data
 */
struct BatchProcessor {
    static constexpr size_t OPTIMAL_BATCH_SIZE = 16;  // Cache-line friendly
    
    std::vector<std::vector<double>> batch_targets;
    std::vector<RangeData> batch_results;
    
    BatchProcessor() {
        batch_targets.reserve(OPTIMAL_BATCH_SIZE);
        batch_results.reserve(OPTIMAL_BATCH_SIZE);
    }
    
    void clear() {
        batch_targets.clear();
        batch_results.clear();
    }
    
    bool isFull() const {
        return batch_targets.size() >= OPTIMAL_BATCH_SIZE;
    }
    
    void addTarget(const std::vector<double>& target_pos) {
        batch_targets.push_back(target_pos);
    }
};

// ============================================================================

int main() {
    // ============================================================================
    // CONFIGURACION DE VERBOSIDAD
    // ============================================================================
    
    // Seleccionar nivel de verbosidad:
    // VerbosityFlags::minimal()  - Solo logs esenciales de progreso y resultados finales
    // VerbosityFlags::standard() - Logs esenciales + resultados numéricos (recomendado)
    // VerbosityFlags::debug()    - Logs completos para debugging (muy verboso)
    
    VerbosityFlags verbosity = VerbosityFlags::minimal(); // Cambiar aquí para ajustar verbosidad
    
    LOG_INFO("Modo de verbosidad: " + std::string(
        verbosity.show_detailed_logs ? "DEBUG" : 
        verbosity.show_numerical_results ? "STANDARD" : "MINIMAL"));
    
    // ============================================================================
    
    // Initialize project directories
    ProjectPaths::ensureDirectoriesExist();
    
    // Inicializar parámetros de radar
    RadarTx strRadarTx = initializeRadarTx();
    RadarRx strRadarRx = initializeRadarRx();
    Target strTarget = initializeTarget();
    System strSystem = initializeSystem();
    Environment strEnvironment = initializeEnvironment();
    
    // Calcular longitud de onda
    strRadarTx.lamb  = strSystem.VelocidadeLuz / strRadarTx.FreqPortadora;
    
    // Cargar datos binarios e interpolar
    std::vector<std::complex<double>> rCompData;
    std::vector<double> grid_xAxis, grid_yAxis, grid_zAxis;
    ProcessedDEM demData;  // Usar estructura compatible
    std::vector<std::vector<double>> Tx, Rx;  // Posiciones transmisor y receptor
    int factorInterp = 4;

    // Actualizar frecuencia de muestreo para receptor
    strRadarRx.fs = strRadarRx.fs * factorInterp;

    try {
        MatrixInfo info;
        std::vector<std::complex<double>> rootData = readBinary<double>(ProjectPaths::getProcPath("rootData_complex.bin"), info);
        LOG_SUCCESS("rootData cargado: " + std::to_string(info.totalElements) + " elementos");
        
        // VERIFICAR: ¿rootData es realmente una matriz 2D?
        LOG_INFO("=== VERIFICACION DE DIMENSIONES rootData ===");
        LOG_INFO("Elementos totales: " + std::to_string(info.totalElements));
        LOG_INFO("Numero de dimensiones: " + std::to_string(info.dimensions.size()));
        LOG_INFO("Dimensiones de la matriz:");
        for (size_t i = 0; i < info.dimensions.size(); ++i) {
            LOG_INFO("  - dimension[" + std::to_string(i) + "]: " + std::to_string(info.dimensions[i]));
        }
        LOG_INFO("Tipo de datos: " + info.dataType);
        LOG_INFO("Es complejo: " + std::string(info.isComplex ? "true" : "false"));
        
        // Verificar si las dimensiones coinciden
        size_t expected_elements = 1;
        for (auto dim : info.dimensions) {
            expected_elements *= dim;
        }
        LOG_INFO("Verificacion: producto de dimensiones = " + std::to_string(expected_elements));
        LOG_INFO("¿Coincide con totalElements? " + std::string(expected_elements == info.totalElements ? "SI" : "NO"));
        
        if (info.dimensions.size() >= 2 && info.dimensions[0] > 1 && info.dimensions[1] > 1) {
            LOG_WARNING("ATENCION: rootData es matriz 2D o multidimensional!");
            LOG_WARNING("Dimensiones: " + std::to_string(info.dimensions[0]) + " x " + std::to_string(info.dimensions[1]));
            LOG_WARNING("La interpolacion y acceso deben considerar estructura 2D");
            
            if (info.dimensions.size() == 2) {
                LOG_INFO("Interpretacion probable: [" + std::to_string(info.dimensions[0]) + " muestras de rango] x [" + 
                         std::to_string(info.dimensions[1]) + " pulsos azimutales]");
            }
        }
        
        // Cargar grids (archivos CSV)
        grid_xAxis = readCSV(ProjectPaths::getProcPath("grid_xAxis.csv"));
        grid_yAxis = readCSV(ProjectPaths::getProcPath("grid_yAxis.csv"));
        grid_zAxis = readCSV(ProjectPaths::getProcPath("grid_zAxis.csv"));

        LOG_SUCCESS("Grids cargados: " + std::to_string(grid_xAxis.size()) + "x" + 
                   std::to_string(grid_yAxis.size()) + "x" + std::to_string(grid_zAxis.size()));
        
        // Obtener dimensiones del grid para cálculos
        size_t nx = grid_xAxis.size();
        size_t ny = grid_yAxis.size(); 
        size_t nz = grid_zAxis.size();
        
        LOG_SUCCESS("Grid dimensions: " + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + " points");
        
        // Cargar DEM en estructura compatible
        demData.X_vec = readCSV(ProjectPaths::getProcPath("DEM_X_vec.csv"));
        demData.Y_vec = readCSV(ProjectPaths::getProcPath("DEM_Y_vec.csv"));
        demData.Z_DEM = readCSVMatrix(ProjectPaths::getProcPath("DEM_Z_matrix.csv"));
        
        LOG_SUCCESS("DEM cargado: X(" + std::to_string(demData.X_vec.size()) + ") Y(" + 
                   std::to_string(demData.Y_vec.size()) + ") Z(" + std::to_string(demData.Z_DEM.size()) + 
                   "x" + std::to_string(demData.Z_DEM.empty() ? 0 : demData.Z_DEM[0].size()) + ")");
        
        // Cargar posiciones Tx y Rx usando función reutilizable
        Tx = readCSVMatrix(ProjectPaths::getProcPath("Tx_positions.csv"));
        Rx = readCSVMatrix(ProjectPaths::getProcPath("Rx_positions.csv"));
        
        LOG_SUCCESS("Posiciones cargadas: Tx(" + std::to_string(Tx.size()) + " posiciones) Rx(" + 
                   std::to_string(Rx.size()) + " posiciones)");
        
        // Preparar DEM para procesamiento usando función compatible
        ProcessedDEM processedDEM = prepareDEMData(demData);
        
        // Variables para procesamiento
        int ratioInterp = 4; // Upsampling ratio for interpolation
        
        // Verificar que tenemos matriz 2D con dimensiones correctas
        if (info.dimensions.size() != 2) {
            LOG_ERROR("Error: rootData debe ser matriz 2D");
            return 1;
        }
        
        size_t numRangeSamples = info.dimensions[0];  // 256 muestras de rango
        size_t numAzimuthPulses = info.dimensions[1]; // 6689 pulsos azimutales
        size_t newRangeSamples = numRangeSamples * factorInterp; // 256 * 4 = 1024
        
        // Interpolar rootData 2D: solo en dimensión fast-time (rango)
        // rootData es matriz 2D: 256 (rango) x 6689 (azimut)
        // Solo interpolamos la dimensión de rango: 256 -> 1024
        LOG_INFO("Datos originales rootData: " + std::to_string(rootData.size()) + " elementos");
        LOG_INFO("Dimensiones originales: " + std::to_string(info.dimensions[0]) + " x " + std::to_string(info.dimensions[1]));
        
        LOG_INFO("Muestras de rango (fast-time): " + std::to_string(numRangeSamples));
        LOG_INFO("Pulsos azimutales (slow-time): " + std::to_string(numAzimuthPulses));
        
        // Crear matriz interpolada 2D: (256*4) x 6689 = 1024 x 6689
        size_t totalInterpolatedElements = newRangeSamples * numAzimuthPulses;
        
        rCompData.resize(totalInterpolatedElements);
        
        LOG_INFO("Interpolando solo en dimensión fast-time...");
        LOG_INFO("Nueva dimensión de rango: " + std::to_string(newRangeSamples));
        LOG_INFO("Preservando dimensión azimutal: " + std::to_string(numAzimuthPulses));
        
        // Interpolar cada pulso azimutal por separado (columna por columna)
        for (size_t azIdx = 0; azIdx < numAzimuthPulses; azIdx++) {
            // Extraer vector de rango para este pulso azimutal
            std::vector<std::complex<double>> rangeVector(numRangeSamples);
            for (size_t rgIdx = 0; rgIdx < numRangeSamples; rgIdx++) {
                size_t originalIdx = rgIdx + azIdx * numRangeSamples; // Row-major order
                rangeVector[rgIdx] = rootData[originalIdx];
            }
            
            // Interpolar este vector de rango
            std::vector<std::complex<double>> interpolatedRange = splineInterpolation(rangeVector, factorInterp);
            
            // Almacenar resultado interpolado en rCompData
            for (size_t rgIdx = 0; rgIdx < newRangeSamples; rgIdx++) {
                size_t newIdx = rgIdx + azIdx * newRangeSamples; // Row-major order
                rCompData[newIdx] = interpolatedRange[rgIdx];
            }
            
            // Progress para interpolación
            if (azIdx % 1000 == 0 || azIdx == numAzimuthPulses - 1) {
                double progress = (double)(azIdx + 1) / numAzimuthPulses * 100.0;
                LOG_PROGRESS("Interpolacion fast-time: " + std::to_string(azIdx + 1) + "/" + 
                           std::to_string(numAzimuthPulses) + " (" + 
                           std::to_string(progress) + "%)");
            }
        }
        
        LOG_SUCCESS("Interpolacion 2D completada: " + std::to_string(rCompData.size()) + " elementos");
        LOG_INFO("Dimensiones interpoladas: " + std::to_string(newRangeSamples) + " x " + std::to_string(numAzimuthPulses));
        LOG_INFO("Factor aplicado solo en fast-time: " + std::to_string(factorInterp) + "x");
        
        // Debug: Mostrar algunos valores de rCompData interpolado
        if (!rCompData.empty()) {
            LOG_INFO("Muestra de rCompData interpolado (primeros 5 elementos):");
            for (size_t i = 0; i < std::min(size_t(5), rCompData.size()); ++i) {
                std::cout << "  [" << i << "] = " << rCompData[i].real() << " + " << rCompData[i].imag() << "i (mag: " << std::abs(rCompData[i]) << ")" << std::endl;
            }
            
            // Verificar valores en diferentes posiciones azimutales para encontrar datos no-cero
            LOG_INFO("Buscando datos no-cero en diferentes posiciones azimutales...");
            bool foundNonZero = false;
            size_t maxSamples = std::min(size_t(10), numAzimuthPulses);
            
            for (size_t azIdx = 0; azIdx < maxSamples && !foundNonZero; azIdx++) {
                size_t baseIdx = azIdx * newRangeSamples;
                for (size_t rgIdx = 340; rgIdx < std::min(baseIdx + newRangeSamples, baseIdx + 350); rgIdx++) {
                    if (baseIdx + rgIdx - baseIdx < rCompData.size()) {
                        size_t actualIdx = baseIdx + (rgIdx - baseIdx);
                        double mag = std::abs(rCompData[actualIdx]);
                        if (mag > 1e-10) { // Threshold para considerar no-cero
                            LOG_INFO("Datos no-cero encontrados en azIdx=" + std::to_string(azIdx) + 
                                   ", rgIdx=" + std::to_string(rgIdx - baseIdx) + 
                                   ", mag=" + std::to_string(mag));
                            foundNonZero = true;
                            break;
                        }
                    }
                }
            }
            
            if (!foundNonZero) {
                LOG_WARNING("No se encontraron datos no-cero en las primeras " + std::to_string(maxSamples) + " posiciones azimutales");
                LOG_WARNING("Los datos pueden estar en otras posiciones azimutales o ser muy pequeños");
            }
        }

        // Calcular distancias usando Fermat para cada punto del grid
        LOG_INFO("Iniciando calculo de rangos Fermat para grid completo...");
        
        // Calcular número total de puntos del grid
        size_t totalGridPoints = nx * ny * nz;
        LOG_INFO("Total de puntos en el grid: " + std::to_string(totalGridPoints));
        
        // Preparar vector target para ser actualizado en cada iteración
        std::vector<double> target_pos(3);
        
        // Variables para almacenar resultados
        std::vector<RangeData> allRanges;
        std::vector<ReflectionPoints> allReflections;
        std::vector<RefractionPoints> allRefractions;
        
        allRanges.reserve(totalGridPoints);
        allReflections.reserve(totalGridPoints);
        allRefractions.reserve(totalGridPoints);
        
        // Crear matriz 3D output para almacenar valores de allPower
        // Dimensiones: nx (X) x ny (Y) x nz (Z)
        // Usar vector anidado para matriz 3D: output[x][y][z]
        std::vector<std::vector<std::vector<std::complex<double>>>> output(
            nx, std::vector<std::vector<std::complex<double>>>(
                ny, std::vector<std::complex<double>>(nz, std::complex<double>(0.0, 0.0))
            )
        );
        
        LOG_INFO("Matriz 3D output creada con dimensiones: " + std::to_string(nx) + " x " + 
                std::to_string(ny) + " x " + std::to_string(nz));
        LOG_INFO("Total de elementos en output: " + std::to_string(nx * ny * nz));
        
        // PARALLEL PROCESSING: Process grid points using OpenMP for real speedup
        LOG_INFO("Iniciando procesamiento paralelo con OpenMP");
        LOG_INFO("Número de threads disponibles: " + std::to_string(omp_get_max_threads()));
        
        // Set OpenMP threads for optimal performance
        omp_set_num_threads(std::min(omp_get_max_threads(), 8)); // Limit to avoid context switching
        
        // Pre-allocate results vectors with proper size
        allRanges.resize(totalGridPoints);
        allReflections.resize(totalGridPoints);
        allRefractions.resize(totalGridPoints);
        
        // Generate flat coordinate list for better cache performance
        std::vector<std::vector<double>> flat_grid_coords;
        generateFlatGridCoordinates(grid_xAxis, grid_yAxis, grid_zAxis, flat_grid_coords);
        
        LOG_INFO("Grid coordinates generadas: " + std::to_string(flat_grid_coords.size()) + " puntos");
        
        // DEBUG: Procesar solo los primeros 8 puntos del grid (comentar las siguientes 3 líneas para procesamiento completo)
        size_t debug_limit = std::min(flat_grid_coords.size(), static_cast<size_t>(8));
        LOG_INFO("DEBUG MODE: Procesando solo los primeros " + std::to_string(debug_limit) + " puntos del grid");
        // Comentar la línea anterior para procesamiento completo
        
        auto parallel_start = std::chrono::high_resolution_clock::now();
        
        // PARALLEL PROCESSING with OpenMP - REAL PERFORMANCE BOOST
        #pragma omp parallel for schedule(dynamic, 1) shared(allRanges, allReflections, allRefractions, flat_grid_coords, Tx, Rx, processedDEM, strEnvironment)
        // for (size_t coord_idx = 0; coord_idx < flat_grid_coords.size(); coord_idx++) {  // Procesamiento completo
        for (size_t coord_idx = 0; coord_idx < debug_limit; coord_idx++) {  // DEBUG: Solo primeros 8 puntos
            
            // Each thread needs its own target structure to avoid race conditions
            Target local_strTarget = initializeTarget();
            
            const auto& target_pos = flat_grid_coords[coord_idx];
            
            // Update local target structure
            updatePointToProc(local_strTarget, target_pos[0], target_pos[1], target_pos[2]);
            
            // Prepare trajectory data for current target position
            TrajectoryMatrix trajectories = prepareTrajectoryData(Tx, Rx, target_pos);
            
            // Log positions for current processing point (MATLAB style) - Thread safe
            #pragma omp critical
            {
                // Get first Tx and Rx positions (assuming they represent the current positions)
                if (!trajectories.Tx.empty() && !trajectories.Rx.empty() && verbosity.show_detailed_logs) {
                    printf("-----------------------------\n");
                    printf("Point %zu - Tx: [%.6f, %.6f, %.6f], Rx: [%.6f, %.6f, %.6f], Pm: [%.6f, %.6f, %.6f]\n",
                           coord_idx + 1,
                           trajectories.Tx[0][0], trajectories.Tx[0][1], trajectories.Tx[0][2],
                           trajectories.Rx[0][0], trajectories.Rx[0][1], trajectories.Rx[0][2],
                           target_pos[0], target_pos[1], target_pos[2]);
                    fflush(stdout); // Force immediate output
                }
            }
            
            // Calculate ranges using Fermat's principle (thread-safe)
            ReflectionPoints reflections;
            RefractionPoints refractions;
            RangeData ranges = calculaSlantRangeFermat(processedDEM, trajectories, strEnvironment, 
                                                     reflections, refractions, false);
            
            // Initialize allPower for this grid point
            std::complex<double> allPower = std::complex<double>(0.0, 0.0);
            
            // Calculate distances using pdist2 equivalent (MATLAB style)
            if (!refractions.P_refrac_ida.empty() && !refractions.P_refrac_volta.empty()) {
                
                // Current point P [1x3] - target_pos
                std::vector<double> P = {target_pos[0], target_pos[1], target_pos[2]};
                
                // R2t = pdist2(P, P_refrac_ida_T) - distances from target to refraction points (ida)
                std::vector<double> R2t;
                R2t.reserve(refractions.P_refrac_ida.size());
                for (const auto& refrac_point : refractions.P_refrac_ida) {
                    double dx = P[0] - refrac_point[0];
                    double dy = P[1] - refrac_point[1]; 
                    double dz = P[2] - refrac_point[2];
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    R2t.push_back(distance);
                }
                
                // R2r = pdist2(P, P_refrac_volta_T) - distances from target to refraction points (volta)
                std::vector<double> R2r;
                R2r.reserve(refractions.P_refrac_volta.size());
                for (const auto& refrac_point : refractions.P_refrac_volta) {
                    double dx = P[0] - refrac_point[0];
                    double dy = P[1] - refrac_point[1];
                    double dz = P[2] - refrac_point[2];
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    R2r.push_back(distance);
                }
                
                // R1t = mean(pdist2(Tx, P_refrac_ida_T), 1) - mean distances from all Tx to refraction points
                std::vector<double> R1t(refractions.P_refrac_ida.size(), 0.0);
                for (size_t i = 0; i < refractions.P_refrac_ida.size(); ++i) {
                    double sum_distances = 0.0;
                    const auto& refrac_point = refractions.P_refrac_ida[i];
                    
                    // Calculate distance from each Tx position to this refraction point
                    for (const auto& tx_pos : trajectories.Tx) {
                        double dx = tx_pos[0] - refrac_point[0];
                        double dy = tx_pos[1] - refrac_point[1];
                        double dz = tx_pos[2] - refrac_point[2];
                        sum_distances += std::sqrt(dx*dx + dy*dy + dz*dz);
                    }
                    R1t[i] = sum_distances / trajectories.Tx.size(); // Mean distance
                }
                
                // R1r = mean(pdist2(Rx, P_refrac_volta_T), 1) - mean distances from all Rx to refraction points  
                std::vector<double> R1r(refractions.P_refrac_volta.size(), 0.0);
                for (size_t i = 0; i < refractions.P_refrac_volta.size(); ++i) {
                    double sum_distances = 0.0;
                    const auto& refrac_point = refractions.P_refrac_volta[i];
                    
                    // Calculate distance from each Rx position to this refraction point
                    for (const auto& rx_pos : trajectories.Rx) {
                        double dx = rx_pos[0] - refrac_point[0];
                        double dy = rx_pos[1] - refrac_point[1];
                        double dz = rx_pos[2] - refrac_point[2];
                        sum_distances += std::sqrt(dx*dx + dy*dy + dz*dz);
                    }
                    R1r[i] = sum_distances / trajectories.Rx.size(); // Mean distance
                }
                
                // Calculate fractional range bin sample and phase compensation (MATLAB equivalent)
                std::vector<double> t, rngBin, phi;
                
                if (!R2t.empty() && !R2r.empty() && !R1t.empty() && !R1r.empty()) {
                    t.reserve(R2t.size());
                    rngBin.reserve(R2t.size());
                    phi.reserve(R2t.size());
                    
                    // Environment refractive indices
                    double n1 = strEnvironment.n1;  // Air refractive index
                    double n2 = strEnvironment.n2;  // Material refractive index
                    double c = strSystem.VelocidadeLuz;  // Speed of light
                    double fs = strRadarRx.fs;  // Sampling frequency
                    double lambda = strRadarTx.lamb;  // Wavelength
                    
                    for (size_t i = 0; i < R2t.size(); ++i) {
                        // Fractional range bin sample: t = (1/c)*(n1*R1t + n2*R2t + n2*R2r + n1*R1r)
                        double t_val = (1.0/c) * (n1*R1t[i] + n2*R2t[i] + n2*R2r[i] + n1*R1r[i]);
                        t.push_back(t_val);
                        
                        // Range bin: rngBin = 1 + t*fs  [MATLAB 1-based indexing, mantener decimales]
                        double rngBin_val = 1.0 + (t_val * fs);
                        rngBin.push_back(rngBin_val);
                        
                        // Phase compensation term: phi = (2*pi/lambda)*(n1*R1t + n2*R2t + n2*R2r + n1*R1r)
                        double phi_val = (2.0 * M_PI / lambda) * (n1*R1t[i] + n2*R2t[i] + n2*R2r[i] + n1*R1r[i]);
                        phi.push_back(phi_val);
                    }
                }
                
                // Back Projection Algorithm - Data Accumulation (MATLAB equivalent)
                // allPower ya está inicializado arriba
                
                if (!rngBin.empty() && !phi.empty() && !rCompData.empty()) {
                    
                    // Para acceso 2D, necesitamos las dimensiones interpoladas
                    size_t interpRangeSamples = newRangeSamples; // 1024 (256 * 4)
                    size_t azimuthPulses = numAzimuthPulses;     // 6689
                    
                    for (size_t idxRng = 0; idxRng < rngBin.size(); ++idxRng) {
                        
                        // MATLAB: idxNext = ceil(rngBin(idxRng)); idxLast = floor(rngBin(idxRng));
                        // Convert from MATLAB 1-based to C++ 0-based indexing
                        double rngBin_val = rngBin[idxRng] - 1.0; // Convert to 0-based
                        
                        int idxLast = static_cast<int>(std::floor(rngBin_val));
                        int idxNext = static_cast<int>(std::ceil(rngBin_val));
                        
                        // Ensure range indices are within bounds
                        if (idxLast >= 0 && idxNext < static_cast<int>(interpRangeSamples)) {
                            
                            // MATLAB: q = rngBin(idxRng) - idxLast;
                            double q = rngBin_val - idxLast;
                            
                            // Para matriz 2D, usamos un índice azimutal específico
                            // En SAR, típicamente corresponde a la posición temporal del pulso
                            // Por simplicidad, usamos una posición azimutal fija para este ejemplo
                            size_t idxAz = std::min(idxRng, numAzimuthPulses - 1); // Mapear a posición azimutal válida
                            
                            // CORRECCION: Acceso correcto a matriz 2D almacenada como 1D
                            // Formato: rCompData[rango + azimut * interpRangeSamples]
                            // donde rango ∈ [0, interpRangeSamples-1] y azimut ∈ [0, numAzimuthPulses-1]
                            
                            // Calcular índices 2D en array 1D (row-major order)
                            size_t idx2D_Last = static_cast<size_t>(idxLast) + idxAz * interpRangeSamples;
                            size_t idx2D_Next = static_cast<size_t>(idxNext) + idxAz * interpRangeSamples;
                            
                            // Verificar que los índices están dentro de los límites
                            bool valid_indices = (idxLast >= 0 && idxNext < static_cast<int>(interpRangeSamples) &&
                                                 idxAz < numAzimuthPulses &&
                                                 idx2D_Last < rCompData.size() && 
                                                 idx2D_Next < rCompData.size());
                            
                            if (valid_indices) {
                                
                                std::complex<double> interpValue;
                                if (idxLast == idxNext) {
                                    // No interpolation needed - acceso directo a matriz 2D
                                    interpValue = rCompData[idx2D_Last];
                                } else {
                                    // Linear interpolation between adjacent range samples
                                    // Ambos valores están en la misma columna azimutal (idxAz)
                                    std::complex<double> val_Last = rCompData[idx2D_Last];
                                    std::complex<double> val_Next = rCompData[idx2D_Next];
                                    interpValue = (1.0 - q) * val_Last + q * val_Next;
                                }
                                
                                // Phase compensation: exp(1i*phi(idxRng))
                                std::complex<double> phaseComp = std::exp(std::complex<double>(0.0, phi[idxRng]));
                                
                                // Accumulate power
                                allPower += interpValue * phaseComp;
                                
                                // Debug: Log some interpolation details for first few iterations
                                if (idxRng < 3 && verbosity.show_debug_info) {
                                    #pragma omp critical
                                    {
                                        printf("    Debug [%zu]: idxAz=%zu, idxLast=%d, idxNext=%d\n", 
                                               idxRng, idxAz, idxLast, idxNext);
                                        printf("    2D indices: idx2D_Last=%zu, idx2D_Next=%zu\n", 
                                               idx2D_Last, idx2D_Next);
                                        printf("    Values: val_Last=%.3e, val_Next=%.3e, interpValue=%.3e\n",
                                               std::abs(rCompData[idx2D_Last]), std::abs(rCompData[idx2D_Next]), std::abs(interpValue));
                                        fflush(stdout);
                                    }
                                }
                            } else {
                                // Debug: Log invalid indices
                                if (idxRng < 3 && verbosity.show_debug_info) {
                                    #pragma omp critical
                                    {
                                        printf("    Warning [%zu]: Invalid indices - idxLast=%d, idxNext=%d, idxAz=%zu\n",
                                               idxRng, idxLast, idxNext, idxAz);
                                        fflush(stdout);
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Log distance calculations and range/phase results (thread-safe)
                #pragma omp critical
                {
                    if (verbosity.show_numerical_results) {
                        printf("---- Distance Calculations (pdist2 equivalent) ----\n");
                        printf("Point %zu - Distance vectors calculated:\n", coord_idx + 1);
                        printf("  R2t size: %zu, R2r size: %zu\n", R2t.size(), R2r.size());
                        printf("  R1t size: %zu, R1r size: %zu\n", R1t.size(), R1r.size());
                        
                        // Show some sample distances
                        if (!R2t.empty() && !R2r.empty() && !R1t.empty() && !R1r.empty() && verbosity.show_detailed_logs) {
                            printf("  Sample distances (first few):\n");
                            for (size_t i = 0; i < std::min(size_t(3), R2t.size()); ++i) {
                                printf("    [%zu] R2t=%.3f, R2r=%.3f, R1t=%.3f, R1r=%.3f\n", 
                                       i, R2t[i], R2r[i], R1t[i], R1r[i]);
                            }
                        }
                        
                        printf("---- Range Bin and Phase Calculations ----\n");
                        printf("Point %zu - Range/Phase vectors calculated:\n", coord_idx + 1);
                        printf("  t size: %zu, rngBin size: %zu, phi size: %zu\n", t.size(), rngBin.size(), phi.size());
                        
                        // Show some sample range/phase calculations
                        if (!t.empty() && !rngBin.empty() && !phi.empty() && verbosity.show_detailed_logs) {
                            printf("  Sample calculations (first few):\n");
                            for (size_t i = 0; i < std::min(size_t(3), t.size()); ++i) {
                                printf("    [%zu] t=%.6e, rngBin=%.3f, phi=%.3f rad\n", 
                                       i, t[i], rngBin[i], phi[i]);
                            }
                        }
                        
                        printf("---- Back Projection Results ----\n");
                        printf("Point %zu - Back Projection Power: Real=%.6f, Imag=%.6f, Magnitude=%.6f\n",
                               coord_idx + 1, allPower.real(), allPower.imag(), std::abs(allPower));
                        
                        // Debug information for range bin access
                        if (!rngBin.empty() && verbosity.show_debug_info) {
                            printf("  Debug info - rCompData size: %zu, rngBin range: [%.3f, %.3f]\n",
                                   rCompData.size(), 
                                   *std::min_element(rngBin.begin(), rngBin.end()),
                                   *std::max_element(rngBin.begin(), rngBin.end()));
                            
                            // Show first few converted indices
                            printf("  Sample converted indices (first 3):\n");
                            for (size_t i = 0; i < std::min(size_t(3), rngBin.size()); ++i) {
                                double rngBin_val = rngBin[i] - 1.0;
                                int idxLast = static_cast<int>(std::floor(rngBin_val));
                                int idxNext = static_cast<int>(std::ceil(rngBin_val));
                                printf("    [%zu] rngBin=%.3f -> rngBin_val=%.3f -> idxLast=%d, idxNext=%d\n",
                                       i, rngBin[i], rngBin_val, idxLast, idxNext);
                            }
                        }
                    }
                    
                    fflush(stdout);
                }
            }
            
            // Log refraction points for specific indices (MATLAB style) - Thread safe
            #pragma omp critical
            {
                if (!refractions.P_refrac_ida.empty() && !refractions.P_refrac_volta.empty() && verbosity.show_refraction_points) {
                    printf("---- Refraction Points ----\n");
                    // Define indices to show (MATLAB: idxToShow = [1,50,100,150,200,250])
                    // Convert from MATLAB 1-based to C++ 0-based: [0,49,99,149,199,249]
                    std::vector<size_t> idxToShow = {0, 49, 99, 149, 199, 249};
                    
                    size_t max_size = std::min(refractions.P_refrac_ida.size(), refractions.P_refrac_volta.size());
                    
                    for (size_t i = 0; i < idxToShow.size(); ++i) {
                        size_t idx = idxToShow[i];
                        if (idx < max_size) {
                            printf("Point %zu - Refracoes %zu: P_refrac_ida: [%.6f, %.6f, %.6f], P_refrac_volta: [%.6f, %.6f, %.6f]\n",
                                   coord_idx + 1, i + 1, // Display as 1-based like MATLAB
                                   refractions.P_refrac_ida[idx][0], refractions.P_refrac_ida[idx][1], refractions.P_refrac_ida[idx][2],
                                   refractions.P_refrac_volta[idx][0], refractions.P_refrac_volta[idx][1], refractions.P_refrac_volta[idx][2]);
                        }
                    }
                    fflush(stdout); // Force immediate output
                }
            }
            
            // Store results (each thread writes to different index - no race condition)
            allRanges[coord_idx] = std::move(ranges);
            allReflections[coord_idx] = std::move(reflections);
            allRefractions[coord_idx] = std::move(refractions);
            
            // Convertir índice lineal coord_idx de vuelta a coordenadas 3D (x,y,z)
            // El orden de generación fue: Z -> X -> Y (Y cambia primero, luego X, luego Z)
            // Por tanto: coord_idx = z*nx*ny + x*ny + y
            size_t y_idx = coord_idx % ny;
            size_t x_idx = (coord_idx / ny) % nx;  
            size_t z_idx = coord_idx / (nx * ny);
            
            // Guardar allPower en la matriz 3D output[x][y][z]
            if (x_idx < nx && y_idx < ny && z_idx < nz) {
                output[x_idx][y_idx][z_idx] = allPower;
                
                // Debug: Log coordinates mapping for first few points
                if (coord_idx < 3 && verbosity.show_coordinate_mapping) {
                    #pragma omp critical
                    {
                        printf("  Coordinate mapping [%zu]: coord_idx=%zu -> x=%zu, y=%zu, z=%zu\n",
                               coord_idx, coord_idx, x_idx, y_idx, z_idx);
                        printf("  Grid position: [%.3f, %.3f, %.3f] -> output[%zu][%zu][%zu] = %.6f\n",
                               target_pos[0], target_pos[1], target_pos[2], 
                               x_idx, y_idx, z_idx, std::abs(allPower));
                        fflush(stdout);
                    }
                }
            }
            
            // Progress logging (thread-safe with reduction)
            #pragma omp critical
            {
                static size_t processed_count = 0;
                processed_count++;
                size_t total_to_process = std::min(flat_grid_coords.size(), static_cast<size_t>(8)); // DEBUG: usar debug_limit
                // size_t total_to_process = totalGridPoints; // Descomentar para procesamiento completo
                if (processed_count % 4 == 0 || processed_count >= total_to_process) { // DEBUG: reporte cada 4 puntos
                    double progress = (double)processed_count / total_to_process * 100.0;
                    LOG_PROGRESS("Parallel processing: " + std::to_string(processed_count) + "/" + 
                               std::to_string(total_to_process) + " (" + 
                               std::to_string(progress) + "%) - Thread: " + 
                               std::to_string(omp_get_thread_num()));
                }
            }
        }
        
        auto parallel_end = std::chrono::high_resolution_clock::now();
        auto parallel_duration = std::chrono::duration_cast<std::chrono::milliseconds>(parallel_end - parallel_start);
        
        size_t points_processed = std::min(flat_grid_coords.size(), static_cast<size_t>(8)); // DEBUG
        // size_t points_processed = totalGridPoints; // Descomentar para procesamiento completo
        
        LOG_PERF("Procesamiento paralelo completado en: " + std::to_string(parallel_duration.count()) + " ms");
        LOG_PERF("Tiempo promedio por punto: " + std::to_string((double)parallel_duration.count() / points_processed) + " ms");
        
        LOG_SUCCESS("Calculo Fermat completado para " + std::to_string(points_processed) + " puntos del grid");
        
        // Generar estadísticas de la matriz 3D output
        if (verbosity.show_numerical_results) {
            LOG_INFO("=== ESTADISTICAS DE LA MATRIZ OUTPUT 3D ===");
            LOG_INFO("Dimensiones: " + std::to_string(nx) + " x " + std::to_string(ny) + " x " + std::to_string(nz));
            
            // Calcular estadísticas básicas
            double min_magnitude = std::numeric_limits<double>::max();
            double max_magnitude = 0.0;
            double sum_magnitude = 0.0;
            size_t non_zero_count = 0;
            
            for (size_t x = 0; x < nx; ++x) {
                for (size_t y = 0; y < ny; ++y) {
                    for (size_t z = 0; z < nz; ++z) {
                        double mag = std::abs(output[x][y][z]);
                        
                        if (mag > 1e-15) { // Considerar no-cero
                            non_zero_count++;
                            min_magnitude = std::min(min_magnitude, mag);
                            sum_magnitude += mag;
                        }
                        max_magnitude = std::max(max_magnitude, mag);
                    }
                }
            }
            
            if (non_zero_count > 0) {
                double avg_magnitude = sum_magnitude / non_zero_count;
                LOG_INFO("Valores no-cero: " + std::to_string(non_zero_count) + "/" + std::to_string(points_processed));
                LOG_INFO("Magnitud mínima: " + std::to_string(min_magnitude));
                LOG_INFO("Magnitud máxima: " + std::to_string(max_magnitude));
                LOG_INFO("Magnitud promedio: " + std::to_string(avg_magnitude));
                
                // Mostrar algunos valores específicos
                if (verbosity.show_detailed_logs) {
                    LOG_INFO("Muestras de valores calculados:");
                    size_t samples_shown = 0;
                    for (size_t x = 0; x < nx && samples_shown < 5; ++x) {
                        for (size_t y = 0; y < ny && samples_shown < 5; ++y) {
                            for (size_t z = 0; z < nz && samples_shown < 5; ++z) {
                                double mag = std::abs(output[x][y][z]);
                                if (mag > 1e-15) {
                                    std::cout << "  output[" << x << "][" << y << "][" << z << "] = " 
                                             << output[x][y][z].real() << " + " << output[x][y][z].imag() 
                                             << "i (mag: " << mag << ")" << std::endl;
                                    samples_shown++;
                                }
                            }
                        }
                    }
                }
            } else {
                LOG_WARNING("No se encontraron valores no-cero en la matriz output");
            }
        } else {
            // Versión simplificada para modo minimal
            size_t non_zero_count = 0;
            for (size_t x = 0; x < nx; ++x) {
                for (size_t y = 0; y < ny; ++y) {
                    for (size_t z = 0; z < nz; ++z) {
                        if (std::abs(output[x][y][z]) > 1e-15) {
                            non_zero_count++;
                        }
                    }
                }
            }
            LOG_SUCCESS("Matriz 3D output completada: " + std::to_string(non_zero_count) + "/" + 
                       std::to_string(points_processed) + " valores no-cero calculados");
        }

        // ============================================================================
        // EXPORTAR MATRIZ 3D OUTPUT
        // ============================================================================
        
        LOG_INFO("=== EXPORTANDO MATRIZ 3D OUTPUT ===");
        
        // Exportar la matriz 3D output a archivo binario
        std::string outputFilename = ProjectPaths::getProcPath("output_3D_matrix.bin");
        export3DComplexMatrix(outputFilename, output, nx, ny, nz, grid_xAxis, grid_yAxis, grid_zAxis);
        
        LOG_SUCCESS("Exportación de matriz 3D completada");
        LOG_INFO("Archivo generado: " + outputFilename);
        
        // ============================================================================

        // Calcular distancias usando Fermat
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error cargando datos: " + std::string(e.what()));
        return 1;
    } 


    return 0;
}