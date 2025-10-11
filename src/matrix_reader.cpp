#include "matrix_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<std::vector<double>> readCSVMatrix(const std::string& filename) {
    std::vector<std::vector<double>> matrix;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir el archivo: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            std::vector<double> row;
            std::stringstream ss(line);
            std::string value;
            
            while (std::getline(ss, value, ',')) {
                if (!value.empty()) {
                    try {
                        row.push_back(std::stod(value));
                    } catch (const std::exception& e) {
                        // Saltar valores que no se pueden convertir
                        continue;
                    }
                }
            }
            
            if (!row.empty()) {
                matrix.push_back(row);
            }
        }
    }
    
    file.close();
    return matrix;
}