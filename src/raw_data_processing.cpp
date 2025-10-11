#include "raw_data_processing.h"
#include "fft_utils.h"
#include "color_logger.h"  // Agregar ColorLogger
#include "project_paths.h" // Centralized paths
#include <iostream>
#include <cmath>
#include <complex>
#include <functional>
#include <algorithm>
#include <chrono>  // Para medición de tiempo
#include <fstream>
#include <sstream>

Environment initializeEnvironment() {
    Environment env;
    
    // Valores por defecto
    env.n1 = 1.0;  // Air refractive index
    env.n2 = 4.0;  // Ground refractive index
    
    // Leer valores desde archivo environment.txt
    std::ifstream file(ProjectPaths::getProcPath("environment.txt"));
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) {
                std::istringstream iss(line);
                std::string key, value;
                
                if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                    if (key == "n1") {
                        env.n1 = std::stod(value);
                    } else if (key == "n2") {
                        env.n2 = std::stod(value);
                    }
                }
            }
        }
        file.close();
        std::cout << "Environment loaded from file: n1=" << env.n1 << ", n2=" << env.n2 << std::endl;
    } else {
        std::cout << "Environment file not found, using defaults: n1=" << env.n1 << ", n2=" << env.n2 << std::endl;
    }
    
    return env;
}

TrajectoryMatrix prepareTrajectoryData(const TrajectoryPoints& trajTx, 
                                      const TrajectoryPoints& trajRx,
                                      const Target& target) {
    TrajectoryMatrix traj;
    
    // Convert TrajectoryPoints to matrix format
    size_t nPoints = trajTx.Px.size();
    traj.Tx.resize(nPoints, std::vector<double>(3));
    traj.Rx.resize(nPoints, std::vector<double>(3));
    
    for (size_t i = 0; i < nPoints; i++) {
        traj.Tx[i][0] = trajTx.Px[i];
        traj.Tx[i][1] = trajTx.Py[i];
        traj.Tx[i][2] = trajTx.Pz[i];
        
        traj.Rx[i][0] = trajRx.Px[i];
        traj.Rx[i][1] = trajRx.Py[i];
        traj.Rx[i][2] = trajRx.Pz[i];
    }
    
    // Set target position
    traj.P.resize(3);
    traj.P[0] = target.pos[0];
    traj.P[1] = target.pos[1];
    traj.P[2] = target.pos[2];
    
    std::cout << "Trajectory data prepared:" << std::endl;
    std::cout << "  Tx positions: " << traj.Tx.size() << std::endl;
    std::cout << "  Rx positions: " << traj.Rx.size() << std::endl;
    std::cout << "  Target: [" << traj.P[0] << ", " << traj.P[1] << ", " << traj.P[2] << "]" << std::endl;
    
    return traj;
}

ProcessedDEM prepareDEMData(const ProcessedDEM& demData) {
    // For this simplified implementation, just return the input data
    // In a more complex version, this could perform coordinate transformations
    std::cout << "DEM data prepared for processing" << std::endl;
    std::cout << "  X vector size: " << demData.X_vec.size() << std::endl;
    std::cout << "  Y vector size: " << demData.Y_vec.size() << std::endl;
    std::cout << "  Z matrix size: " << demData.Z_DEM.size() << "x" 
              << (demData.Z_DEM.empty() ? 0 : demData.Z_DEM[0].size()) << std::endl;
    
    return demData;
}

TrajectoryMatrix prepareTrajectoryData(const std::vector<std::vector<double>>& Tx,
                                      const std::vector<std::vector<double>>& Rx,
                                      const std::vector<double>& target_pos) {
    TrajectoryMatrix traj;
    traj.Tx = Tx;
    traj.Rx = Rx;
    traj.P = target_pos;
    
    std::cout << "Trajectory data prepared:" << std::endl;
    std::cout << "  Tx positions: " << Tx.size() << std::endl;
    std::cout << "  Rx positions: " << Rx.size() << std::endl;
    std::cout << "  Target: [" << target_pos[0] << ", " << target_pos[1] << ", " << target_pos[2] << "]" << std::endl;
    
    return traj;
}

RangeData calculaSlantRangeFermat(const ProcessedDEM& demData,
                                 const TrajectoryMatrix& trajectories,
                                 const Environment& env,
                                 ReflectionPoints& strReflexao,
                                 RefractionPoints& strRefraccoes,
                                 bool bPlotVerbose) {
    
    // [TIMER] Iniciar medicion de tiempo
    auto start_time = std::chrono::high_resolution_clock::now();
    
    LOG_INFO("> Calculating ranges using Fermat's principle...");
    
    size_t nPoints = trajectories.Tx.size();
    RangeData ranges;
    
    // Initialize range vectors
    ranges.r1Tx.resize(nPoints);
    ranges.r2Tx.resize(nPoints);
    ranges.r1Rx.resize(nPoints);
    ranges.r2Rx.resize(nPoints);
    
    // Initialize refraction point structures
    strRefraccoes.P_refrac_ida.resize(nPoints, std::vector<double>(3));
    strRefraccoes.P_refrac_volta.resize(nPoints, std::vector<double>(3));
    strReflexao.P_reflection.resize(nPoints, std::vector<double>(3));
    
    // Target position
    double Px = trajectories.P[0];
    double Py = trajectories.P[1]; 
    double Pz = trajectories.P[2];
    
    if (bPlotVerbose) {
        std::cout << "Fermat calculation with verbose output enabled" << std::endl;
        std::cout << "Target position: [" << Px << ", " << Py << ", " << Pz << "]" << std::endl;
        std::cout << "Refractive indices: n1=" << env.n1 << ", n2=" << env.n2 << std::endl;
    }
    
    // Process each trajectory point
    auto start_fermat_loop = std::chrono::high_resolution_clock::now();
    
    for (size_t i = 0; i < nPoints; i++) {
        // Tx position
        double TxX = trajectories.Tx[i][0];
        double TxY = trajectories.Tx[i][1];
        double TxZ = trajectories.Tx[i][2];
        
        // Rx position  
        double RxX = trajectories.Rx[i][0];
        double RxY = trajectories.Rx[i][1];
        double RxZ = trajectories.Rx[i][2];
        
        // Calculate refraction point for Tx->Target path (ida)
        // This implements Snell's law at the interface to find optimal refraction point
        std::vector<double> P_refrac_ida = calculateRefractionPoint(
            TxX, TxY, TxZ,    // Transmitter position
            Px, Py, Pz,      // Target position
            demData,          // DEM surface data
            env.n1, env.n2,  // Refractive indices
            bPlotVerbose      // Debug output flag
        );
        
        // Calculate refraction point for Target->Rx path (volta)
        std::vector<double> P_refrac_volta = calculateRefractionPoint(
            Px, Py, Pz,      // Target position
            RxX, RxY, RxZ,   // Receiver position
            demData,          // DEM surface data
            env.n2, env.n1,  // Refractive indices (reversed for return path)
            bPlotVerbose      // Debug output flag
        );
        
        // Store refraction points
        strRefraccoes.P_refrac_ida[i] = P_refrac_ida;
        strRefraccoes.P_refrac_volta[i] = P_refrac_volta;
        
        // Calculate reflection point (simplified - at DEM surface)
        std::vector<double> P_reflection = calculateReflectionPoint(
            TxX, TxY, TxZ,
            RxX, RxY, RxZ,
            demData
        );
        strReflexao.P_reflection[i] = P_reflection;
        
        // Calculate ranges
        // r1Tx: Distance from Tx to refraction point (ida)
        ranges.r1Tx[i] = sqrt(
            pow(TxX - P_refrac_ida[0], 2) +
            pow(TxY - P_refrac_ida[1], 2) +
            pow(TxZ - P_refrac_ida[2], 2)
        );
        
        // r2Tx: Distance from refraction point (ida) to target
        ranges.r2Tx[i] = sqrt(
            pow(Px - P_refrac_ida[0], 2) +
            pow(Py - P_refrac_ida[1], 2) +
            pow(Pz - P_refrac_ida[2], 2)
        );
        
        // r1Rx: Distance from Rx to refraction point (volta)
        ranges.r1Rx[i] = sqrt(
            pow(RxX - P_refrac_volta[0], 2) +
            pow(RxY - P_refrac_volta[1], 2) +
            pow(RxZ - P_refrac_volta[2], 2)
        );
        
        // r2Rx: Distance from target to refraction point (volta)
        ranges.r2Rx[i] = sqrt(
            pow(Px - P_refrac_volta[0], 2) +
            pow(Py - P_refrac_volta[1], 2) +
            pow(Pz - P_refrac_volta[2], 2)
        );
        
        if (bPlotVerbose) {
            std::cout << "Point " << i << ": Tx[" << TxX << "," << TxY << "," << TxZ << "] -> "
                      << "Refrac_ida[" << P_refrac_ida[0] << "," << P_refrac_ida[1] << "," << P_refrac_ida[2] << "] -> "
                      << "Target[" << Px << "," << Py << "," << Pz << "]" << std::endl;

            std::cout << "Point " << i << ": Target[" << Px << "," << Py << "," << Pz << "] -> "
                      << "Refrac_volta[" << P_refrac_volta[0] << "," << P_refrac_volta[1] << "," << P_refrac_volta[2] << "] -> "
                      << "Rx[" << RxX << "," << RxY << "," << RxZ << "]" << std::endl;

            std::cout << "  Ranges (m): r1Tx=" << ranges.r1Tx[i] 
                      << ", r2Tx=" << ranges.r2Tx[i]
                      << ", r1Rx=" << ranges.r1Rx[i]
                      << ", r2Rx=" << ranges.r2Rx[i] << std::endl;

            std::cout << " -----------------------------------------------------------------------------------" << std::endl;
        }
    }
    
    auto end_fermat_loop = std::chrono::high_resolution_clock::now();
    auto fermat_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_fermat_loop - start_fermat_loop);
    
    // [TIMER] Finalizar medicion de tiempo y mostrar resultados
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    LOG_SUCCESS("+ Fermat range calculations completed for " + std::to_string(nPoints) + " positions");
    LOG_TIMER("T Total execution time: " + std::to_string(duration.count()) + " ms");
    LOG_INFO("i Fermat loop time: " + std::to_string(fermat_duration.count()) + " ms (" + std::to_string(fermat_duration.count() * 100.0 / duration.count()) + "% of total)");
    LOG_PROGRESS("* Average time per point: " + std::to_string(fermat_duration.count() * 1000.0 / nPoints) + " us/point");
              
    return ranges;
}

std::vector<std::complex<double>> createCHIRP(const RadarTx& radarTx) {
    std::cout << "Creating CHIRP signal..." << std::endl;
    
    // Extract parameters exactly as in MATLAB version
    double f0 = radarTx.FreqMenor;     // f0 = strRadar.FreqMenor;
    double f1 = radarTx.FreqMaior;     // f1 = strRadar.FreqMayor;
    double T = radarTx.DuracaoPulso;   // T  = strRadar.DuracaoPulso;
    double fs = radarTx.fs;            // fs = strRadar.fs;
    
    // Calculate chirp rate: k = (f1 - f0)/T;
    double k = (f1 - f0) / T;
    
    // Generate time vector: t = 0:1/fs:T;
    // Number of samples including endpoint
    int N = static_cast<int>(T * fs) + 1;
    std::vector<std::complex<double>> chirp(N);
    
    // Generate chirp: chirp = exp(1j*2*pi*(f0*t + k/2*t.^2));
    for (int n = 0; n < N; n++) {
        double t = n / fs;  // Time at sample n
        
        // Calculate phase: 2*pi*(f0*t + k/2*t^2)
        double phase = 2.0 * M_PI * (f0 * t + (k / 2.0) * t * t);
        
        // Create complex exponential: exp(1j*phase)
        chirp[n] = std::complex<double>(cos(phase), sin(phase));
    }
    
    double B = f1 - f0;  // Bandwidth for display
    LOG_SUCCESS("+ CHIRP signal created: " + std::to_string(N) + " samples, B=" + std::to_string(B/1e6) + " MHz");
    std::cout << "  f0=" << f0/1e6 << " MHz, f1=" << f1/1e6 << " MHz, T=" << T*1e6 << " μs" << std::endl;
    
    return chirp;
}

RawDataResult processRawData(const TrajectoryMatrix& trajectories,
                            const RangeData& ranges,
                            const RadarTx& radarTx,
                            const RadarRx& radarRx,
                            const Target& target,
                            const System& system,
                            const Environment& env) {
    std::cout << "Processing raw radar data..." << std::endl;
    
    RawDataResult result;
    
    double c = system.VelocidadeLuz;
    size_t nPoints = trajectories.Tx.size();
    
    // Create CHIRP signal: chirp = createCHRIP(strRadarTx);
    result.chirp = createCHIRP(radarTx);
    
    // Initialize data matrices
    result.auxData.resize(system.IndiceMaximo, std::vector<std::complex<double>>(nPoints, 0.0));
    result.rawData.resize(system.IndiceMaximo, std::vector<std::complex<double>>(nPoints, 0.0));
    result.rootData.resize(system.IndiceMaximo, std::vector<std::complex<double>>(nPoints, 0.0));
    
    std::cout << "Data matrices initialized: " << system.IndiceMaximo << " x " << nPoints << std::endl;
    
    // Step 1: Generate auxData - equivalent to MATLAB lines 130-135
    // phi = -(2*pi/lamb)*(n1*r1Tx + n2*r2Tx + n1*r1Rx + n2*r2Rx);
    // auxData = zeros(strSystem.IndiceMaximo,length(PxT));
    // IND = sub2ind(size(auxData),rngBin,1:length(PxT));
    // auxData(IND) = strTarget.rcs.*exp(1i*phi);
    
    for (size_t i = 0; i < nPoints; i++) {
        // Calculate total optical path
        double totalRange = env.n1 * ranges.r1Tx[i] + env.n2 * ranges.r2Tx[i] + 
                           env.n1 * ranges.r1Rx[i] + env.n2 * ranges.r2Rx[i];
        double t = totalRange / c;
        
        // Calculate range bin: rngBin = 1 + round(t * fs)
        int rngBin = 1 + static_cast<int>(round(t * radarRx.fs));
        
        // Alert if total range seems excessive (much larger than direct path)
        if (i > 0) {
            double directTxTarget = sqrt(pow(trajectories.Tx[i][0] - trajectories.P[0], 2) + 
                                       pow(trajectories.Tx[i][1] - trajectories.P[1], 2) + 
                                       pow(trajectories.Tx[i][2] - trajectories.P[2], 2));
            double directTargetRx = sqrt(pow(trajectories.P[0] - trajectories.Rx[i][0], 2) + 
                                       pow(trajectories.P[1] - trajectories.Rx[i][1], 2) + 
                                       pow(trajectories.P[2] - trajectories.Rx[i][2], 2));
            double directPath = directTxTarget + directTargetRx;
            
            if (totalRange > directPath * 1.5) {  // If Fermat path is >50% longer than direct
                LOG_WARNING("! Point " + std::to_string(i) + ": Excessive path length detected");
                LOG_WARNING("  Total optical: " + std::to_string(totalRange) + "m, Direct: " + std::to_string(directPath) + "m");
                LOG_WARNING("  r1Tx=" + std::to_string(ranges.r1Tx[i]) + ", r2Tx=" + std::to_string(ranges.r2Tx[i]) + 
                           ", r1Rx=" + std::to_string(ranges.r1Rx[i]) + ", r2Rx=" + std::to_string(ranges.r2Rx[i]));
            }
        }
        
        // Calculate phase: phi = -(2*pi/lamb) * totalRange
        double phi = -(2.0 * M_PI / radarTx.lamb) * totalRange;
        
        // Place target response: auxData(IND) = strTarget.rcs.*exp(1i*phi);
        if (rngBin >= 0 && rngBin < system.IndiceMaximo) {
            result.auxData[rngBin][i] = target.rcs * std::complex<double>(cos(phi), sin(phi));
        }

        // Show t, rngBin and phi for first few points for verification (one line)
        // if (nPoints <= 10 || i < 10) {
            LOG_DEBUG("? Point " + std::to_string(i) + ": t=" + std::to_string(t*1e9) + " ns, rngBin=" + std::to_string(rngBin) + ", phi=" + std::to_string(phi) + " rad");
        // }
    }
    
    LOG_SUCCESS("+ Target responses placed in auxData matrix");
    
    // Step 2: Create reference signal - equivalent to MATLAB lines 136-142
    // pulsoTx = zeros(1,strSystem.IndiceMaximo);
    // pulsoTx(1:length(chirp)) = chirp;
    // K = fix(length(chirp)/2);
    // reference = circshift([pulsoTx],K).';
    
    std::vector<std::complex<double>> pulsoTx(system.IndiceMaximo, 0.0);
    
    // Fill with chirp: pulsoTx(1:length(chirp)) = chirp;
    size_t chirpLength = result.chirp.size();
    for (size_t i = 0; i < std::min(chirpLength, (size_t)system.IndiceMaximo); i++) {
        pulsoTx[i] = result.chirp[i];
    }
    
    // FFTShift pattern: K = fix(length(chirp)/2); [segunda_mitad, zeros, primera_mitad]
    int K = chirpLength / 2;
    std::vector<std::complex<double>> reference(system.IndiceMaximo, 0.0);
    
    // Implement fftshift pattern like in our successful test:
    // Place second half of chirp at the beginning of the padded signal
    int half_chirp = chirpLength / 2;
    for (int i = 0; i < (chirpLength - half_chirp); ++i) {
        reference[i] = pulsoTx[half_chirp + i];
    }
    
    // Place first half of chirp at the end of the padded signal
    int start_pos_first_half = system.IndiceMaximo - half_chirp;
    for (int i = 0; i < half_chirp; ++i) {
        reference[start_pos_first_half + i] = pulsoTx[i];
    }
    
    LOG_SUCCESS("+ Reference signal created with fftshift pattern [" + std::to_string(chirpLength - half_chirp) + " second_half + " + 
                std::to_string(start_pos_first_half - (chirpLength - half_chirp)) + " zeros + " + 
                std::to_string(half_chirp) + " first_half] = " + std::to_string(system.IndiceMaximo) + " total samples");
    
    // Step 3: Raw data processing - equivalent to MATLAB line 143
    // rawData = ifft(fft(auxData).*fft(reference));
    
    std::cout << "Processing raw data with FFT correlation..." << std::endl;
    
    // For each column (azimuth position)
    for (size_t col = 0; col < nPoints; col++) {
        // Extract column from auxData
        std::vector<std::complex<double>> auxCol(system.IndiceMaximo);
        for (int row = 0; row < system.IndiceMaximo; row++) {
            auxCol[row] = result.auxData[row][col];
        }
        
        // Perform convolution: ifft(fft(auxCol) .* fft(reference)) - WITHOUT conjugate
        // Note: This matches MATLAB line 143: rawData = ifft(fft(auxData).*fft(reference));
        std::vector<std::complex<double>> rawCol = FFTProcessor::fftConvolution(auxCol, reference);
        
        // Store result
        for (int row = 0; row < system.IndiceMaximo; row++) {
            result.rawData[row][col] = rawCol[row];
        }
    }
    
    LOG_SUCCESS("Raw data correlation completed");
    
    // Step 4: Pulse compression - equivalent to MATLAB line 146
    // rootData = ifft(fft(rawData).*conj(fft(reference)));
    
    std::cout << "Applying pulse compression..." << std::endl;
    
    // For each column (azimuth position)
    for (size_t col = 0; col < nPoints; col++) {
        // Extract column from rawData
        std::vector<std::complex<double>> rawCol(system.IndiceMaximo);
        for (int row = 0; row < system.IndiceMaximo; row++) {
            rawCol[row] = result.rawData[row][col];
        }
        
        // Perform compression: ifft(fft(rawCol) .* conj(fft(reference))) - WITH conjugate
        // Note: This matches MATLAB line 146: rootData = ifft(fft(rawData).*conj(fft(reference)));
        std::vector<std::complex<double>> rootCol = FFTProcessor::fftCorrelation(rawCol, reference);
        
        // Store result
        for (int row = 0; row < system.IndiceMaximo; row++) {
            result.rootData[row][col] = rootCol[row];
        }
    }
    
    LOG_SUCCESS("● Pulse compression completed");
    LOG_SUCCESS("● Raw data processing completed - equivalent to MATLAB lines 138-146");
    
    return result;
}

// Helper function to calculate refraction point using Fermat's principle with optimization
std::vector<double> calculateRefractionPoint(double x1, double y1, double z1,   // Start point
                                            double x2, double y2, double z2,   // End point
                                            const ProcessedDEM& demData,       // DEM surface data
                                            double n1, double n2,              // Refractive indices
                                            bool bPlotVerbose) {               // Debug output flag
    
    std::vector<double> refractionPoint(3);
    
    // Always use the full Fermat optimization algorithm
    // Even without valid DEM data, we'll use a flat surface at z=0 for the optimization
    if (bPlotVerbose) {
        if (demData.X_vec.empty() || demData.Y_vec.empty() || demData.Z_DEM.empty()) {
            LOG_DEBUG("? No DEM data available, using flat surface at z=0 for Fermat optimization");
        } else {
            LOG_DEBUG("? Using DEM surface data for Fermat optimization");
        }
        LOG_DEBUG("? Executing full Fermat algorithm:");
    }
    
    // Simplified optimization using golden section search in 2D
    // This replaces the full Nelder-Mead but still finds a good minimum
    
    // Initial search bounds - for flat surface, search more extensively initially
    double midX = (x1 + x2) / 2.0;
    double midY = (y1 + y2) / 2.0;
    
    // For flat surface, expand search radius to cover a larger area initially
    // This ensures we find the global minimum, not a local one
    double directDistance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    double searchRadius = std::max(directDistance, 500.0);  // Start with larger search area (min 500m)
    
    if (bPlotVerbose) {
        LOG_DEBUG("? Fermat search initialization:");
        LOG_DEBUG("  Start point: [" + std::to_string(x1) + ", " + std::to_string(y1) + ", " + std::to_string(z1) + "]");
        LOG_DEBUG("  End point: [" + std::to_string(x2) + ", " + std::to_string(y2) + ", " + std::to_string(z2) + "]");
        LOG_DEBUG("  Search center: [" + std::to_string(midX) + ", " + std::to_string(midY) + "]");
        LOG_DEBUG("  Direct distance: " + std::to_string(directDistance) + " m");
        LOG_DEBUG("  Search radius: " + std::to_string(searchRadius) + " m");
        LOG_DEBUG("  Refractive indices: n1=" + std::to_string(n1) + ", n2=" + std::to_string(n2));
    }
    
    double bestX = midX, bestY = midY;
    double minOpticalPath = 1e9;
    
    // Grid search with refinement - start broad, then narrow down
    int numLevels = 5;  // More refinement levels for better accuracy
    double currentRadius = searchRadius;
    
    if (bPlotVerbose) {
        LOG_DEBUG("  Starting grid search with " + std::to_string(numLevels) + " refinement levels");
    }
    
    for (int level = 0; level < numLevels; level++) {
        int gridSize = (level == 0) ? 15 : 11;  // Use finer grid for first level
        double minX = bestX - currentRadius;
        double maxX = bestX + currentRadius;
        double minY = bestY - currentRadius;
        double maxY = bestY + currentRadius;
        
        if (bPlotVerbose) {
            LOG_DEBUG("  Level " + std::to_string(level) + ": grid=" + std::to_string(gridSize) + "x" + std::to_string(gridSize) + 
                     ", radius=" + std::to_string(currentRadius) + " m");
            LOG_DEBUG("    Search bounds: X[" + std::to_string(minX) + ", " + std::to_string(maxX) + "], Y[" + 
                     std::to_string(minY) + ", " + std::to_string(maxY) + "]");
        }
        
        double levelBestOptical = 1e9;
        
        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                double px = minX + i * (maxX - minX) / (gridSize - 1);
                double py = minY + j * (maxY - minY) / (gridSize - 1);
                
                // Get surface elevation at point (px, py)
                double pz = interpolateDEMElevation(px, py, demData);
                
                // Calculate optical path: n1*d1 + n2*d2 (Fermat's principle)
                double d1 = sqrt(pow(px-x1, 2) + pow(py-y1, 2) + pow(pz-z1, 2));
                double d2 = sqrt(pow(x2-px, 2) + pow(y2-py, 2) + pow(z2-pz, 2));
                double opticalPath = n1 * d1 + n2 * d2;
                
                if (opticalPath < minOpticalPath) {
                    minOpticalPath = opticalPath;
                    bestX = px;
                    bestY = py;
                    levelBestOptical = opticalPath;
                }
            }
        }
        
        if (bPlotVerbose) {
            LOG_DEBUG("    Best at level " + std::to_string(level) + ": [" + std::to_string(bestX) + ", " + 
                     std::to_string(bestY) + "], optical_path=" + std::to_string(levelBestOptical) + " m");
        }
        
        // Refine search radius for next level - more aggressive refinement
        currentRadius /= (level == 0) ? 5.0 : 3.0;
    }
    
    // Set final refraction point
    refractionPoint[0] = bestX;
    refractionPoint[1] = bestY;
    refractionPoint[2] = interpolateDEMElevation(bestX, bestY, demData);
    
    if (bPlotVerbose) {
        double finalD1 = sqrt(pow(bestX-x1, 2) + pow(bestY-y1, 2) + pow(refractionPoint[2]-z1, 2));
        double finalD2 = sqrt(pow(x2-bestX, 2) + pow(y2-bestY, 2) + pow(z2-refractionPoint[2], 2));
        double finalOpticalPath = n1 * finalD1 + n2 * finalD2;
        
        LOG_DEBUG("? Fermat optimization result:");
        LOG_DEBUG("  Final refraction point: [" + std::to_string(bestX) + ", " + std::to_string(bestY) + ", " + std::to_string(refractionPoint[2]) + "]");
        LOG_DEBUG("  Segment distances: d1=" + std::to_string(finalD1) + " m, d2=" + std::to_string(finalD2) + " m");
        LOG_DEBUG("  Final optical path: " + std::to_string(finalOpticalPath) + " m");
        LOG_DEBUG("  Path efficiency: " + std::to_string(finalOpticalPath/directDistance * 100.0) + "% of direct distance");
    }
    
    return refractionPoint;
}

// Helper function to calculate reflection point (simplified)
std::vector<double> calculateReflectionPoint(double txX, double txY, double txZ,
                                           double rxX, double rxY, double rxZ,
                                           const ProcessedDEM& demData) {
    
    std::vector<double> reflectionPoint(3);
    
    // Simplified reflection point calculation
    // Find the midpoint projected onto the DEM surface
    double midX = (txX + rxX) / 2.0;
    double midY = (txY + rxY) / 2.0;
    
    // Find closest DEM point to midpoint
    double minDist = 1e9;
    double bestZ = 0.0;
    
    if (!demData.X_vec.empty() && !demData.Y_vec.empty() && !demData.Z_DEM.empty()) {
        for (size_t i = 0; i < demData.X_vec.size(); i++) {
            for (size_t j = 0; j < demData.Y_vec.size(); j++) {
                if (j < demData.Z_DEM.size() && i < demData.Z_DEM[j].size()) {
                    double dx = demData.X_vec[i] - midX;
                    double dy = demData.Y_vec[j] - midY;
                    double dist = sqrt(dx*dx + dy*dy);
                    
                    if (dist < minDist) {
                        minDist = dist;
                        bestZ = demData.Z_DEM[j][i];
                    }
                }
            }
        }
    }
    
    reflectionPoint[0] = midX;
    reflectionPoint[1] = midY;
    reflectionPoint[2] = bestZ;
    
    return reflectionPoint;
}

// Helper function to interpolate DEM elevation at any (x,y) point
double interpolateDEMElevation(double x, double y, const ProcessedDEM& demData) {
    if (demData.X_vec.empty() || demData.Y_vec.empty() || demData.Z_DEM.empty()) {
        return 0.0;  // Default ground level
    }
    
    // Find surrounding grid points for bilinear interpolation
    size_t i0 = 0, i1 = 0;
    size_t j0 = 0, j1 = 0;
    
    // Find X indices
    for (size_t i = 0; i < demData.X_vec.size() - 1; i++) {
        if (x >= demData.X_vec[i] && x <= demData.X_vec[i + 1]) {
            i0 = i;
            i1 = i + 1;
            break;
        }
    }
    
    // Find Y indices
    for (size_t j = 0; j < demData.Y_vec.size() - 1; j++) {
        if (y >= demData.Y_vec[j] && y <= demData.Y_vec[j + 1]) {
            j0 = j;
            j1 = j + 1;
            break;
        }
    }
    
    // Clamp to valid range
    i0 = std::min(i0, demData.X_vec.size() - 1);
    i1 = std::min(i1, demData.X_vec.size() - 1);
    j0 = std::min(j0, demData.Y_vec.size() - 1);
    j1 = std::min(j1, demData.Y_vec.size() - 1);
    
    // Get corner elevations (check bounds)
    double z00 = (j0 < demData.Z_DEM.size() && i0 < demData.Z_DEM[j0].size()) ? demData.Z_DEM[j0][i0] : 0.0;
    double z01 = (j0 < demData.Z_DEM.size() && i1 < demData.Z_DEM[j0].size()) ? demData.Z_DEM[j0][i1] : 0.0;
    double z10 = (j1 < demData.Z_DEM.size() && i0 < demData.Z_DEM[j1].size()) ? demData.Z_DEM[j1][i0] : 0.0;
    double z11 = (j1 < demData.Z_DEM.size() && i1 < demData.Z_DEM[j1].size()) ? demData.Z_DEM[j1][i1] : 0.0;
    
    // Bilinear interpolation
    if (i0 == i1 && j0 == j1) {
        return z00;  // Point exact match
    }
    
    double dx = (i1 < demData.X_vec.size() && i0 < demData.X_vec.size()) ? 
                (demData.X_vec[i1] - demData.X_vec[i0]) : 1.0;
    double dy = (j1 < demData.Y_vec.size() && j0 < demData.Y_vec.size()) ? 
                (demData.Y_vec[j1] - demData.Y_vec[j0]) : 1.0;
    
    if (dx == 0) dx = 1.0;
    if (dy == 0) dy = 1.0;
    
    double wx = (i1 < demData.X_vec.size()) ? (x - demData.X_vec[i0]) / dx : 0.0;
    double wy = (j1 < demData.Y_vec.size()) ? (y - demData.Y_vec[j0]) / dy : 0.0;
    
    wx = std::max(0.0, std::min(1.0, wx));
    wy = std::max(0.0, std::min(1.0, wy));
    
    // Bilinear interpolation formula
    double z = z00 * (1 - wx) * (1 - wy) +
               z01 * wx * (1 - wy) +
               z10 * (1 - wx) * wy +
               z11 * wx * wy;
    
    return z;
}

// Helper function to perform FFT-based correlation using efficient Cooley-Tukey algorithm
// Equivalent to MATLAB: ifft(fft(signal) .* fft(reference))
std::vector<std::complex<double>> correlateWithReference(
    const std::vector<std::complex<double>>& signal,
    const std::vector<std::complex<double>>& reference) {
        
    // Use our efficient FFT implementation
    std::vector<std::complex<double>> result = FFTProcessor::fftCorrelation(signal, reference);
    
    return result;
}