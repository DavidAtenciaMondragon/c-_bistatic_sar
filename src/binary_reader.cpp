#include "binary_reader.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <algorithm>

template<typename T>
std::vector<std::complex<T>> readBinary(const std::string& filename, MatrixInfo& info) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir el archivo: " + filename);
    }
    
    try {
        // 2. Leer el Encabezado
        
        // 2.1. Número de dimensiones (uint32)
        uint32_t numDimensions;
        file.read(reinterpret_cast<char*>(&numDimensions), sizeof(numDimensions));
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudo leer el número de dimensiones");
        }
        
        // 2.2. Bandera de complejidad (uint32): 0=Real, 1=Complejo
        uint32_t isComplexFlag;
        file.read(reinterpret_cast<char*>(&isComplexFlag), sizeof(isComplexFlag));
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudo leer la bandera de complejidad");
        }
        info.isComplex = (isComplexFlag == 1);
        
        // 2.3. Longitud del tipo de dato (uint32) y la cadena (char[16] fijo)
        uint32_t dataTypeLen;
        file.read(reinterpret_cast<char*>(&dataTypeLen), sizeof(dataTypeLen));
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudo leer la longitud del tipo de dato");
        }
        
        char fixedDataTypeStr[16];
        file.read(fixedDataTypeStr, 16);
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudo leer la cadena del tipo de dato");
        }
        
        // Extraer el nombre real del tipo de dato, eliminando los NULLs de relleno
        info.dataType = std::string(fixedDataTypeStr, std::min(static_cast<size_t>(dataTypeLen), static_cast<size_t>(16)));
        
        // 2.4. Dimensiones (uint64[N])
        info.dimensions.resize(numDimensions);
        file.read(reinterpret_cast<char*>(info.dimensions.data()), numDimensions * sizeof(uint64_t));
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudieron leer las dimensiones");
        }
        
        // 3. Calcular el número total de valores a leer
        info.totalElements = 1;
        for (uint64_t dim : info.dimensions) {
            info.totalElements *= dim;
        }
        
        size_t numValuesToRead = info.totalElements * (info.isComplex ? 2 : 1);
        
        // 4. Leer los datos de la matriz (BLOB)
        std::vector<T> matrixValues(numValuesToRead);
        file.read(reinterpret_cast<char*>(matrixValues.data()), numValuesToRead * sizeof(T));
        if (!file.good()) {
            throw std::runtime_error("Error de lectura: No se pudieron leer todos los datos de la matriz");
        }
        
        // 5. Reconstruir la matriz
        std::vector<std::complex<T>> result;
        result.reserve(info.totalElements);
        
        if (info.isComplex) {
            // Los datos se leyeron como [R1, I1, R2, I2, ...].
            // Convertir a números complejos
            for (size_t i = 0; i < matrixValues.size(); i += 2) {
                T real = matrixValues[i];
                T imag = matrixValues[i + 1];
                result.emplace_back(real, imag);
            }
        } else {
            // Datos reales, convertir a complejos con parte imaginaria cero
            for (size_t i = 0; i < matrixValues.size(); i++) {
                result.emplace_back(matrixValues[i], T(0));
            }
        }
        
        // Log de éxito (simplificado)
        // std::cout << "✅ Matriz leída exitosamente: " << filename << std::endl;
        
        return result;
        
    } catch (const std::exception& e) {
        file.close();
        throw;
    }
}

// Especialización explícita para double
template std::vector<std::complex<double>> readBinary<double>(const std::string& filename, MatrixInfo& info);

// Función de conveniencia para leer solo datos reales
template<typename T>
std::vector<T> readRealBinary(const std::string& filename, MatrixInfo& info) {
    auto complexData = readBinary<T>(filename, info);
    std::vector<T> realData;
    realData.reserve(complexData.size());
    
    for (const auto& val : complexData) {
        realData.push_back(val.real());
    }
    
    return realData;
}

// Especialización explícita para double
template std::vector<double> readRealBinary<double>(const std::string& filename, MatrixInfo& info);