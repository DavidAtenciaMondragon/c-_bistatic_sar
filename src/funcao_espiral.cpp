#include "funcao_espiral.h"
#include <cmath>
#include <vector>

TrajectoryPoints funcao_espiral(double N, double R_minor, double R_major, 
                               double H_high, double H_low, double Vt, double PRF) {
    
    // Start with R_major and end with R_minor (as requested)
    double R_start = R_major;  // Always start with major radius
    double R_end = R_minor;    // Always end with minor radius
    
    // Calculate derived parameters
    double Vr = Vt * (std::log(R_start) - std::log(R_end)) / (N * 2 * M_PI);  // Radius reduction rate [m/s]
    double T = (R_start - R_end) / Vr;                                          // Total time [s]
    double Vd = (H_high - H_low) / T;                                           // Descent velocity [m/s]
    
    // Create time vector
    std::vector<double> t;
    double dt = 1.0 / PRF;  // Time step [s]
    double current_time = 0.0;
    
    while (current_time <= T) {
        t.push_back(current_time);
        current_time += dt;
    }
    
    // Initialize result structure
    TrajectoryPoints result;
    result.Px.reserve(t.size());
    result.Py.reserve(t.size());
    result.Pz.reserve(t.size());
    
    // Calculate trajectory points
    for (size_t i = 0; i < t.size(); i++) {
        double time = t[i];
        
        // Calculate angle [rad]
        double theta = (Vt / Vr) * (std::log(R_start) - std::log(R_start - Vr * time));
        
        // Calculate radius [m]
        double R = R_start - Vr * time;
        
        // Calculate coordinates
        result.Px.push_back(R * std::cos(theta));
        result.Py.push_back(R * std::sin(theta));
        result.Pz.push_back(H_high - Vd * time);
    }
    
    return result;
}