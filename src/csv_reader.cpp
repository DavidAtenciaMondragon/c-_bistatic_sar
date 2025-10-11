#include "csv_reader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<double> readCSV(const std::string& filename) {
    std::vector<double> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("No se pudo abrir archivo CSV: " + filename);
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                if (!value.empty()) {
                    data.push_back(std::stod(value));
                }
            }
        }
    }
    return data;
}