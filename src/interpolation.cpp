#include "interpolation.h"
#include <cmath>

std::vector<std::complex<double>> splineInterpolation(const std::vector<std::complex<double>>& data, int factor) {
    size_t n = data.size();
    size_t newSize = n * factor;
    std::vector<std::complex<double>> result(newSize);
    
    // Para datos complejos, interpolar parte real e imaginaria por separado
    std::vector<double> realPart(n), imagPart(n);
    for (size_t i = 0; i < n; i++) {
        realPart[i] = data[i].real();
        imagPart[i] = data[i].imag();
    }
    
    // Interpolación spline simple usando interpolación lineal suavizada
    for (size_t i = 0; i < newSize; i++) {
        double pos = static_cast<double>(i) / factor;
        size_t idx = static_cast<size_t>(pos);
        double frac = pos - idx;
        
        if (idx >= n - 1) {
            idx = n - 2;
            frac = 1.0;
        }
        
        // Interpolación cúbica simple (Hermite)
        double t = frac;
        double t2 = t * t;
        double t3 = t2 * t;
        
        // Coeficientes de Hermite
        double h1 = 2*t3 - 3*t2 + 1;
        double h2 = -2*t3 + 3*t2;
        double h3 = t3 - 2*t2 + t;
        double h4 = t3 - t2;
        
        // Calcular derivadas (diferencias finitas)
        double realDerivPrev = (idx > 0) ? (realPart[idx] - realPart[idx-1]) : 0;
        double realDerivNext = (idx < n-2) ? (realPart[idx+2] - realPart[idx+1]) : 0;
        double imagDerivPrev = (idx > 0) ? (imagPart[idx] - imagPart[idx-1]) : 0;
        double imagDerivNext = (idx < n-2) ? (imagPart[idx+2] - imagPart[idx+1]) : 0;
        
        // Interpolación spline
        double realInterp = h1 * realPart[idx] + h2 * realPart[idx+1] + h3 * realDerivPrev + h4 * realDerivNext;
        double imagInterp = h1 * imagPart[idx] + h2 * imagPart[idx+1] + h3 * imagDerivPrev + h4 * imagDerivNext;
        
        result[i] = std::complex<double>(realInterp, imagInterp);
    }
    
    return result;
}