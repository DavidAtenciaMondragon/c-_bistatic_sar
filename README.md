# Radar Trajectory C++ Project

Este proyecto traduce el script MATLAB `GS_fermat_test_script.m` a C++, incluyendo las funciones de trayectoria espiral y carga de archivos DEM.

## Dependencias

### GDAL (Geospatial Data Abstraction Library)

Para compilar este proyecto, necesitas instalar GDAL:

#### Windows:
1. Descargar GDAL desde: https://www.gisinternals.com/
2. O usar vcpkg:
   ```powershell
   vcpkg install gdal
   ```

#### Ubuntu/Debian:
```bash
sudo apt-get install libgdal-dev
```

#### macOS:
```bash
brew install gdal
```

## Compilación

1. Crear directorio de construcción:
   ```bash
   mkdir build
   cd build
   ```

2. Configurar con CMake:
   ```bash
   cmake ..
   ```

3. Compilar:
   ```bash
   cmake --build .
   ```

## Uso

1. Asegúrate de que el archivo DEM esté en la ruta correcta: `assets/DEM_1x1km_Res30m_Lat-3_6160_Lon-80_4552.tif`

2. Ejecutar:
   ```bash
   ./bin/radar_trajectory
   ```

## Estructura del proyecto

```
├── include/
│   ├── radar_params.h      # Estructuras de parámetros del radar
│   ├── funcao_espiral.h    # Función de trayectoria espiral
│   └── load_dem.h          # Función de carga de DEM
├── src/
│   ├── funcao_espiral.cpp  # Implementación de trayectoria espiral
│   └── load_dem.cpp        # Implementación de carga de DEM
├── main.cpp                # Programa principal
└── CMakeLists.txt          # Configuración de compilación
```

## Funciones traducidas

- ✅ `funcao_espiral()` - Generación de trayectoria espiral
- ✅ `loadDEM()` - Carga de archivos GeoTIFF con conversión UTM
- ✅ `plotEscenario()` - Generación de datos de visualización (equivalente al plot MATLAB)
- ✅ Carga de parámetros del radar
- ✅ Cálculo de longitud de onda

## Visualización

El programa genera archivos CSV con los datos del escenario que puedes visualizar:

### Archivos generados en `plot_data/`:
- `transmitter_trajectory.csv` - Trayectoria del transmisor (línea negra en MATLAB)
- `receiver_trajectory.csv` - Trayectoria del receptor (línea roja en MATLAB)
- `target_position.csv` - Posición del target (marcador negro en MATLAB)
- `dem_surface.csv` - Superficie DEM (con transparencia en MATLAB)
- `plot_scenario.py` - Script de Python para visualización automática

### Para visualizar:

**Opción 1: Python (automático)**
```bash
cd build/plot_data
python plot_scenario.py
```

**Opción 2: MATLAB**
```matlab
% Cargar datos
tx = readtable('plot_data/transmitter_trajectory.csv');
rx = readtable('plot_data/receiver_trajectory.csv');
target = readtable('plot_data/target_position.csv');
dem = readtable('plot_data/dem_surface.csv');

% Plot equivalente al original
figure
hold on
plot3(tx.X, tx.Y, tx.Z, 'k')
plot3(rx.X, rx.Y, rx.Z, 'r')
plot3(target.X, target.Y, target.Z, 'o', 'MarkerFaceColor', 'k')
% Para el DEM necesitas reshape según dem_grid_info.txt
legend("Tx","Rx","Target")
```

## Notas

- Los parámetros están codificados en las funciones `initialize*()`. 
- Para usar valores reales, reemplaza estos valores con datos de tus archivos JSON.
- El archivo DEM debe estar en formato GeoTIFF (.tif)