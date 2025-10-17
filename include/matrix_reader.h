#ifndef MATRIX_READER_H
#define MATRIX_READER_H

#include <vector>
#include <string>

/**
 * @brief Lee un archivo CSV que contiene una matriz de valores double
 * @param filename Ruta al archivo CSV
 * @return Matriz 2D (vector de vectores) con los datos del archivo
 * @throws std::runtime_error si no se puede abrir el archivo
 */
std::vector<std::vector<double>> readCSVMatrix(const std::string& filename);

#endif // MATRIX_READER_H