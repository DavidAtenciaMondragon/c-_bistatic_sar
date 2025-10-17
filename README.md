# Sistema SAR Bistático - Documentación Completa

## 🎯 Descripción General

Sistema completo de procesamiento de radar SAR (Synthetic Aperture Radar) bistático implementado en C++ con capacidades avanzadas de:

- ✅ **Generación de trayectorias espirales** para plataformas Tx/Rx
- ✅ **Procesamiento de datos DEM** (Digital Elevation Model)
- ✅ **Cálculos de Fermat** para propagación con refracción
- ✅ **Algoritmo de back projection** para formación de imágenes SAR
- ✅ **Exportación de matrices 3D complejas** en formato binario optimizado
- ✅ **Visualización Python** para análisis de resultados
- ✅ **Sistema de verbosidad configurable** para logs
- ✅ **Paralelización OpenMP** para alto rendimiento

## 📁 Estructura del Proyecto

```
core/
├── 📋 CMakeLists.txt              # Configuración de build
├── 🔨 compile.bat                 # Script de compilación Windows
├── 🧹 clean_build.bat             # Script de limpieza
├── 📖 README.md                   # Este archivo
├── 📊 PROJECT_STRUCTURE.md        # Documentación de organización
│
├── 📁 include/                    # Headers (.h)
│   ├── binary_reader.h           # Lectura archivos binarios
│   ├── color_logger.h            # Sistema de logs con colores
│   ├── csv_reader.h              # Lectura archivos CSV
│   ├── data_export.h             # ⭐ Exportación datos (incluye 3D)
│   ├── interpolation.h           # Algoritmos interpolación
│   ├── matrix_reader.h           # Lectura matrices
│   ├── radar_params.h            # Estructuras parámetros radar
│   ├── raw_data_processing.h     # Procesamiento datos raw
│   └── ... (otros headers)
│
├── 📁 src/                       # Implementaciones (.cpp)
│   ├── data_export.cpp           # ⭐ Exportación 3D matrices
│   ├── raw_data_processing.cpp   # Algoritmos SAR core
│   ├── interpolation.cpp         # Fast-time interpolation
│   └── ... (otras implementaciones)
│
├── 📁 Scripts principales
│   ├── 🚀 proc_test_script.cpp    # ⭐ Back projection + export 3D
│   ├── 🔬 gs_test_script.cpp      # Generación datos + export
│   └── 📊 example_colored_logs.cpp # Ejemplo sistema logs
│
├── 📁 Visualización
│   ├── 🎨 visualize_binary_data.py # Visualizador datos binarios
│   └── 🎯 visualize_3D_matrix.py   # ⭐ Visualizador matriz 3D SAR
│
└── 📁 build/                     # Archivos de build
    ├── bin/                      # Ejecutables
    └── output/                   # Datos procesados
        └── proc/                 # Archivos exportados
            ├── ⭐ output_3D_matrix.bin    # Matriz 3D SAR
            ├── output_3D_matrix_grid_*.csv # Coordenadas grid
            ├── rootData_complex.bin      # Datos radar complejos
            ├── Tx_positions.csv          # Posiciones transmisor
            ├── Rx_positions.csv          # Posiciones receptor
            └── environment.txt           # Parámetros ambiente
```

## 🚀 Compilación y Ejecución

### Requisitos
- **CMake 3.10+**
- **Ninja** (build system)
- **GCC/G++ 13.2+** con soporte C++17
- **OpenMP** para paralelización
- **Python 3.x** + numpy + matplotlib (para visualización)

### Compilación Rápida
```bash
# Opción 1: Script automático (Windows)
.\compile.bat

# Opción 2: Comandos manuales
cmake -S . -B build -G Ninja
ninja -C build
```

### Ejecución

#### 1. Procesamiento SAR Completo (⭐ Recomendado)
```bash
# Ejecutar back projection con exportación 3D
.\build\bin\proc_test_script.exe

# Salida esperada:
# - Matriz 3D: build/output/proc/output_3D_matrix.bin
# - Coordenadas: build/output/proc/output_3D_matrix_grid_*.csv
# - Logs con estadísticas de procesamiento
```

#### 2. Generación de Datos Base
```bash
# Generar datos de entrada para procesamiento
.\build\bin\gs_test_script.exe

# Salida esperada:
# - rootData_complex.bin (datos radar)
# - Tx_positions.csv, Rx_positions.csv (trayectorias)
# - DEM_*.csv (modelo elevación)
```

#### 3. Visualización de Resultados
```bash
# Analizar matriz 3D SAR generada
python visualize_3D_matrix.py

# O especificar archivo específico:
python visualize_3D_matrix.py build/output/proc/output_3D_matrix.bin
```

## ⚙️ Configuración del Sistema

### Sistema de Verbosidad

El script principal (`proc_test_script.cpp`) incluye niveles configurables de verbosidad:

```cpp
// En proc_test_script.cpp línea ~184:
VerbosityFlags verbosity = VerbosityFlags::standard(); // ⬅️ Cambiar aquí

// Opciones disponibles:
VerbosityFlags::minimal()   // Solo logs esenciales + resultado final
VerbosityFlags::standard()  // Logs + resultados numéricos (recomendado)
VerbosityFlags::debug()     // Logs completos para debugging (muy verboso)
```

### Modo Debug vs Producción

```cpp
// En proc_test_script.cpp líneas ~434-436:
size_t debug_limit = std::min(flat_grid_coords.size(), static_cast<size_t>(8));
LOG_INFO("DEBUG MODE: Procesando solo los primeros " + std::to_string(debug_limit) + " puntos del grid");

// Para procesamiento completo (35,301 puntos):
// 1. Comentar las líneas de debug_limit
// 2. Descomentar líneas de procesamiento completo
// 3. Recompilar
```

## 📊 Formato de Datos

### Matriz 3D SAR (output_3D_matrix.bin)

**Formato binario optimizado:**

```
Header:
├── numDimensions (uint32)     # 3 para matriz 3D
├── isComplex (uint32)         # 1 = datos complejos
├── dataTypeLen (uint32)       # Longitud string tipo
├── dataType (char[16])        # "double" (fijo 16 bytes)
├── dimensions (uint64[3])     # [nx, ny, nz]
└── data (double[])            # Real-Imaginario intercalado

Datos:
└── Complex values orden C: matrix[x][y][z] = real + i*imag
```

**Dimensiones típicas:**
- **41 × 41 × 21** = 35,301 elementos complejos
- **Tamaño archivo:** ~565 KB por matriz
- **Coordenadas reales:** X[19.2, 20.8]m, Y[-0.8, 0.8]m, Z[-5.8, -4.2]m

### Coordenadas del Grid

Archivos CSV asociados con coordenadas espaciales:
- `output_3D_matrix_grid_x.csv` - Coordenadas X (metros)
- `output_3D_matrix_grid_y.csv` - Coordenadas Y (metros)  
- `output_3D_matrix_grid_z.csv` - Coordenadas Z (metros)

## 🎯 Flujo de Procesamiento SAR

### 1. Inicialización
```
📡 Parámetros Radar → 🛰️ Trayectorias → 🗺️ DEM → 📐 Grid 3D
```

### 2. Procesamiento de Datos
```
📊 rootData (256×6689) → 🔄 Interpolación Fast-time → 📈 rCompData (1024×6689)
```

### 3. Back Projection Algorithm
```
🎯 Para cada punto del grid:
├── 📏 Calcular distancias Fermat (pdist2 equivalent)
├── ⏱️ Calcular range bins y compensación fase
├── 🔄 Interpolación bilinear en datos 2D
├── 📡 Aplicar compensación fase
└── 📊 Acumular power en matriz 3D
```

### 4. Exportación y Análisis
```
💾 Matriz 3D → 📁 Archivo binario → 🎨 Visualización Python
```

## 📈 Resultados Típicos

### Estadísticas de Procesamiento (Modo Debug)
```
✅ Grid: 41×41×21 = 35,301 puntos totales
✅ Procesados: 8 puntos (modo debug)
✅ Valores no-cero: 8/8 (100% en debug)
✅ Rango magnitudes: [91.97, 771.92]
✅ Tiempo promedio: ~1.25 segundos/punto
✅ Paralelización: 8 threads OpenMP
```

### Visualización 3D
- **Histograma de magnitudes**
- **Proyecciones 2D** (XY, XZ, YZ)
- **Cortes transversales**
- **Scatter plot 3D** de puntos significativos
- **Localización de máximos**

## 🔧 Desarrollo y Extensión

### Agregar Nueva Funcionalidad
1. **Header:** Agregar declaración en `/include/nombre.h`
2. **Implementación:** Crear `/src/nombre.cpp`
3. **CMakeLists:** Agregar a las dependencias apropiadas
4. **Testing:** Usar modo debug para validación

### Optimización de Performance
- **Paralelización:** Ajustar `omp_set_num_threads()`
- **Memoria:** Optimizar tamaño de batch en `BatchProcessor`
- **I/O:** Usar SSD para archivos temporales grandes

### Debugging
```cpp
// Activar logs detallados
VerbosityFlags verbosity = VerbosityFlags::debug();

// Activar banderas específicas
verbosity.show_debug_info = true;          // Info interpolación
verbosity.show_refraction_points = true;   // Puntos refracción  
verbosity.show_coordinate_mapping = true;  // Mapeo coordenadas
```

## 📋 Troubleshooting

### Errores Comunes

**Error de Compilación:**
```bash
# Limpiar build cache
rm -r build/CMakeFiles
cmake -S . -B build -G Ninja
ninja -C build
```

**Error "matrix_reader.h not found":**
```
✅ Solucionado: Archivo movido a /include/ (ver PROJECT_STRUCTURE.md)
```

**Datos no-cero faltantes:**
```cpp
// Verificar interpolación y acceso 2D
verbosity.show_debug_info = true;
```

**Performance lenta:**
```cpp
// Reducir threads si hay context switching
omp_set_num_threads(4); // En lugar de 8+
```

## 🌟 Características Avanzadas

### ⭐ **Exportación 3D Optimizada**
- Formato binario compacto
- Header autodescriptivo
- Compatibilidad multiplataforma
- Metadatos de coordenadas integrados

### ⭐ **Visualización Integrada**
- Análisis estadístico completo
- Múltiples tipos de proyección
- Localización automática de máximos
- Exportación de gráficos HD

### ⭐ **Sistema de Logs Inteligente**
- Colores configurables por tipo
- Niveles de verbosidad granulares
- Thread-safe para paralelización
- Separadores estilo MATLAB

## 📞 Soporte

Para preguntas técnicas o reportar problemas:
1. Revisar logs con `VerbosityFlags::debug()`
2. Verificar estructura de archivos con `PROJECT_STRUCTURE.md`
3. Probar con datasets reducidos en modo debug
4. Consultar visualizaciones para validar resultados

---

**🎯 Sistema SAR Bistático** - Implementación completa C++ con exportación 3D y visualización Python integrada.