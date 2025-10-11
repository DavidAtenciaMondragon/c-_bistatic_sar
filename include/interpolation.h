#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <complex>

// Función para interpolación spline cúbica de datos complejos
std::vector<std::complex<double>> splineInterpolation(const std::vector<std::complex<double>>& data, int factor);

#endif // INTERPOLATION_H