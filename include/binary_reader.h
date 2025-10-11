#ifndef BINARY_READER_H
#define BINARY_READER_H

#include <vector>
#include <complex>
#include <string>
#include <cstdint>

// Estructura para almacenar información de la matriz leída
struct MatrixInfo {
    std::vector<uint64_t> dimensions;
    bool isComplex;
    std::string dataType;
    size_t totalElements;
};

// Función principal para leer matrices binarias
template<typename T>
std::vector<std::complex<T>> readComplexBinary(const std::string& filename, MatrixInfo& info);

template<typename T>
std::vector<T> readRealBinary(const std::string& filename, MatrixInfo& info);

// Función genérica que determina automáticamente si es real o complejo
template<typename T>
std::vector<std::complex<T>> readBinary(const std::string& filename, MatrixInfo& info);

#endif // BINARY_READER_H