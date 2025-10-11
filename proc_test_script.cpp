#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include "radar_params.h"
#include "radar_init.h"
#include "trajectory_generator.h"
#include "binary_reader.h"
#include "interpolation.h" 
#include "csv_reader.h"
#include "src/matrix_reader.h"
#include "raw_data_processing.h"
#include "color_logger.h"
#include "project_paths.h" // Centralized paths
#include <omp.h>  // Para paralelización OpenMP

// Color logger para output mejorado
#define LOG_INFO(msg) std::cout << "[INFO] " << msg << std::endl
#define LOG_ERROR(msg) std::cerr << "[ERROR] " << msg << std::endl  
#define LOG_PROGRESS(msg) std::cout << "[PROGRESS] " << msg << std::endl
#define LOG_PERF(msg) std::cout << "[PERFORMANCE] " << msg << std::endl

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
    for (const auto& z : grid_z) {
        for (const auto& y : grid_y) {
            for (const auto& x : grid_x) {
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
        
        // Cargar grids (archivos CSV)
        grid_xAxis = readCSV(ProjectPaths::getProcPath("grid_xAxis.csv"));
        grid_yAxis = readCSV(ProjectPaths::getProcPath("grid_yAxis.csv"));
        grid_zAxis = readCSV(ProjectPaths::getProcPath("grid_zAxis.csv"));
        LOG_SUCCESS("Grids cargados: " + std::to_string(grid_xAxis.size()) + "x" + 
                   std::to_string(grid_yAxis.size()) + "x" + std::to_string(grid_zAxis.size()));
        
        // Crear meshgrid equivalente a [X,Y,Z] = meshgrid(grid_xAxis,grid_yAxis,grid_zAxis)
        size_t nx = grid_xAxis.size();
        size_t ny = grid_yAxis.size(); 
        size_t nz = grid_zAxis.size();
        
        // Crear matrices 3D para X, Y, Z
        std::vector<std::vector<std::vector<double>>> X(ny, std::vector<std::vector<double>>(nx, std::vector<double>(nz)));
        std::vector<std::vector<std::vector<double>>> Y(ny, std::vector<std::vector<double>>(nx, std::vector<double>(nz)));
        std::vector<std::vector<std::vector<double>>> Z(ny, std::vector<std::vector<double>>(nx, std::vector<double>(nz)));
        
        // Llenar las matrices meshgrid
        for (size_t i = 0; i < ny; i++) {        // grid_yAxis (primer índice en MATLAB)
            for (size_t j = 0; j < nx; j++) {    // grid_xAxis (segundo índice en MATLAB)
                for (size_t k = 0; k < nz; k++) { // grid_zAxis (tercer índice en MATLAB)
                    X[i][j][k] = grid_xAxis[j];  // X varía con j (columnas)
                    Y[i][j][k] = grid_yAxis[i];  // Y varía con i (filas)
                    Z[i][j][k] = grid_zAxis[k];  // Z varía con k (profundidad)
                }
            }
        }
        
        LOG_SUCCESS("Meshgrid creado: " + std::to_string(ny) + "x" + std::to_string(nx) + "x" + std::to_string(nz) + " matrices X,Y,Z");
        
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
        
        // Interpolar rootData usando interpolación spline cúbica
        rCompData = splineInterpolation(rootData, factorInterp);
        
        LOG_SUCCESS("Interpolacion completada: " + std::to_string(rCompData.size()) + " elementos");

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
        
        auto parallel_start = std::chrono::high_resolution_clock::now();
        
        // PARALLEL PROCESSING with OpenMP - REAL PERFORMANCE BOOST
        #pragma omp parallel for schedule(dynamic, 1) shared(allRanges, allReflections, allRefractions, flat_grid_coords, Tx, Rx, processedDEM, strEnvironment)
        for (size_t coord_idx = 0; coord_idx < flat_grid_coords.size(); coord_idx++) {
            
            // Each thread needs its own target structure to avoid race conditions
            Target local_strTarget = initializeTarget();
            
            const auto& target_pos = flat_grid_coords[coord_idx];
            
            // Update local target structure
            updatePointToProc(local_strTarget, target_pos[0], target_pos[1], target_pos[2]);
            
            // Prepare trajectory data for current target position
            TrajectoryMatrix trajectories = prepareTrajectoryData(Tx, Rx, target_pos);
            
            // Calculate ranges using Fermat's principle (thread-safe)
            ReflectionPoints reflections;
            RefractionPoints refractions;
            RangeData ranges = calculaSlantRangeFermat(processedDEM, trajectories, strEnvironment, 
                                                     reflections, refractions, false);
            
            // Store results (each thread writes to different index - no race condition)
            allRanges[coord_idx] = std::move(ranges);
            allReflections[coord_idx] = std::move(reflections);
            allRefractions[coord_idx] = std::move(refractions);
            
            // Progress logging (thread-safe with reduction)
            #pragma omp critical
            {
                static size_t processed_count = 0;
                processed_count++;
                if (processed_count % 100 == 0 || processed_count >= totalGridPoints) {
                    double progress = (double)processed_count / totalGridPoints * 100.0;
                    LOG_PROGRESS("Parallel processing: " + std::to_string(processed_count) + "/" + 
                               std::to_string(totalGridPoints) + " (" + 
                               std::to_string(progress) + "%) - Thread: " + 
                               std::to_string(omp_get_thread_num()));
                }
            }
        }
        
        auto parallel_end = std::chrono::high_resolution_clock::now();
        auto parallel_duration = std::chrono::duration_cast<std::chrono::milliseconds>(parallel_end - parallel_start);
        
        LOG_PERF("Procesamiento paralelo completado en: " + std::to_string(parallel_duration.count()) + " ms");
        LOG_PERF("Tiempo promedio por punto: " + std::to_string((double)parallel_duration.count() / totalGridPoints) + " ms");
        
        LOG_SUCCESS("Calculo Fermat completado para " + std::to_string(totalGridPoints) + " puntos del grid");

        // Calcular distancias usando Fermat
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error cargando datos: " + std::string(e.what()));
        return 1;
    } 


    return 0;
}