# Estructura del Proyecto - Bistatic SAR Processing

## Organización de Archivos

### `/include/` - Headers (Archivos .h)
Contiene todas las declaraciones de funciones, clases y estructuras:
- `binary_reader.h` - Lectura de archivos binarios
- `color_logger.h` - Sistema de logs con colores
- `crear_grid.h` - Generación de grids
- `csv_reader.h` - Lectura de archivos CSV
- `data_export.h` - Exportación de datos
- `fft_utils.h` - Utilidades FFT
- `funcao_espiral.h` - Funciones espirales
- `interpolation.h` - Algoritmos de interpolación
- `load_dem.h` - Carga de datos DEM
- `matrix_reader.h` - **[CORREGIDO]** Lectura de matrices (movido desde /src)
- `plot_escenario.h` - Visualización de escenarios
- `project_paths.h` - Gestión de rutas centralizadas
- `radar_init.h` - Inicialización de parámetros radar
- `radar_params.h` - Estructuras de parámetros radar
- `raw_data_processing.h` - Procesamiento de datos raw
- `trajectory_generator.h` - Generación de trayectorias

### `/src/` - Implementaciones (Archivos .cpp)
Contiene todas las implementaciones correspondientes a los headers:
- `binary_reader.cpp`
- `crear_grid.cpp`
- `csv_reader.cpp`
- `data_export.cpp`
- `fft_utils.cpp`
- `funcao_espiral.cpp`
- `interpolation.cpp`
- `load_dem.cpp`, `load_dem_simple.cpp`
- `matrix_reader.cpp`
- `plot_escenario.cpp`
- `project_paths.cpp`
- `radar_init.cpp`
- `raw_data_processing.cpp`
- `trajectory_generator.cpp`

### Archivos Principales
- `proc_test_script.cpp` - Programa principal con algoritmo de back projection
- `gs_test_script.cpp` - Script de prueba alternativo
- `example_colored_logs.cpp` - Ejemplo de uso del sistema de logs

### Sistema de Build
- `CMakeLists.txt` - Configuración de CMake
- `compile.bat` - Script de compilación para Windows
- `clean_build.bat` - Script de limpieza

## Convenciones de Organización

### ✅ **Correcto**
- Headers (.h) en `/include/`
- Implementaciones (.cpp) en `/src/`
- Un header por cada archivo .cpp correspondiente
- Includes relativos desde la raíz del proyecto

### ❌ **Incorrecto (corregido)**
- ~~Headers en `/src/` (como estaba `matrix_reader.h`)~~
- ~~Includes con rutas absolutas desde `/src/`~~

## Sistema de Verbosidad

El código principal incluye un sistema de flags de verbosidad configurable:

```cpp
// Niveles disponibles:
VerbosityFlags::minimal()   // Solo logs esenciales
VerbosityFlags::standard()  // Logs + resultados numéricos (recomendado)
VerbosityFlags::debug()     // Logs completos para debugging
```

## Notas de Desarrollo

- Todos los headers deben estar en `/include/` para consistencia
- Los includes deben usar rutas relativas desde la raíz del proyecto
- El sistema de build (CMake) está configurado para buscar headers en `/include/`
- La compilación con ninja se realiza desde `/build/`