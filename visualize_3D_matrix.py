#!/usr/bin/env python3
"""
Visualizador de Matriz 3D SAR - Output Processor
=================================================

Script independiente para leer y visualizar matrices 3D complejas exportadas
desde el sistema de procesamiento SAR bistático.

Archivo: visualize_3D_matrix.py
Autor: Sistema SAR Bistático C++
Fecha: 2025-10-17

Uso:
    python visualize_3D_matrix.py [archivo_3D.bin]

Funcionalidades:
- Lee matrices 3D complejas en formato binario
- Visualiza distribuciones de magnitud
- Genera mapas de calor 2D por capas Z
- Exporta estadísticas detalladas
- Localiza máximos locales en el espacio 3D
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import struct
import sys
import os

def read_3D_complex_matrix(filename):
    """
    Lee matriz 3D compleja desde archivo binario con formato específico.
    
    Formato del archivo:
    1. numDimensions (uint32) - Número de dimensiones (esperado: 3)
    2. isComplex (uint32) - Flag de complejidad (1=complejo, 0=real)
    3. dataTypeLen (uint32) - Longitud del string de tipo de datos
    4. dataType (char[16]) - Tipo de datos fijo de 16 caracteres
    5. dimensions (uint64[3]) - Dimensiones X, Y, Z
    6. data (double[]) - Datos en formato R-I intercalado
    
    Returns:
        tuple: (matrix_3d, dimensions, info_dict)
    """
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Archivo no encontrado: {filename}")
    
    print(f"📁 Leyendo matriz 3D desde: {filename}")
    
    with open(filename, 'rb') as f:
        # 1. Leer número de dimensiones
        num_dimensions = struct.unpack('<I', f.read(4))[0]
        print(f"   Dimensiones del archivo: {num_dimensions}")
        
        if num_dimensions != 3:
            raise ValueError(f"Se esperaban 3 dimensiones, encontradas: {num_dimensions}")
        
        # 2. Leer flag de complejidad
        is_complex = struct.unpack('<I', f.read(4))[0]
        print(f"   Es complejo: {'Sí' if is_complex else 'No'}")
        
        if not is_complex:
            raise ValueError("Se esperaba matriz compleja")
        
        # 3. Leer longitud del tipo de datos
        data_type_len = struct.unpack('<I', f.read(4))[0]
        print(f"   Longitud tipo datos: {data_type_len}")
        
        # 4. Leer tipo de datos (char[16] fijo)
        data_type_bytes = f.read(16)
        data_type = data_type_bytes.decode('utf-8').rstrip('\x00')
        print(f"   Tipo de datos: '{data_type}'")
        
        # 5. Leer dimensiones (3 uint64)
        dims = struct.unpack('<QQQ', f.read(24))  # 3 × 8 bytes
        nx, ny, nz = dims
        total_elements = nx * ny * nz
        
        print(f"   Dimensiones matriz: {nx} × {ny} × {nz}")
        print(f"   Elementos totales: {total_elements}")
        
        # 6. Leer datos complejos (R-I intercalado)
        print("   📊 Leyendo datos complejos...")
        
        # Cada elemento complejo = 2 doubles (16 bytes)
        expected_bytes = total_elements * 16
        remaining_bytes = os.path.getsize(filename) - f.tell()
        
        print(f"   Bytes esperados: {expected_bytes}")
        print(f"   Bytes disponibles: {remaining_bytes}")
        
        if remaining_bytes < expected_bytes:
            raise ValueError(f"Archivo truncado: faltan {expected_bytes - remaining_bytes} bytes")
        
        # Leer todos los datos de una vez
        raw_data = f.read(expected_bytes)
        
        # Convertir a array de doubles (formato little-endian)
        doubles_array = np.frombuffer(raw_data, dtype='<f8')  # f8 = double
        
        # Verificar cantidad de datos
        expected_doubles = total_elements * 2  # Real + Imaginario por elemento
        actual_doubles = len(doubles_array)
        
        print(f"   Doubles esperados: {expected_doubles}")
        print(f"   Doubles leídos: {actual_doubles}")
        
        if actual_doubles != expected_doubles:
            raise ValueError(f"Error en datos: esperados {expected_doubles}, leídos {actual_doubles}")
        
        # Convertir a números complejos
        # Formato: [R1, I1, R2, I2, R3, I3, ...]
        real_parts = doubles_array[0::2]  # Elementos pares (índices 0, 2, 4, ...)
        imag_parts = doubles_array[1::2]  # Elementos impares (índices 1, 3, 5, ...)
        
        complex_array = real_parts + 1j * imag_parts
        
        # Reshape a matriz 3D (x varía primero, luego y, luego z)
        matrix_3d = complex_array.reshape((nx, ny, nz), order='C')
        
        print(f"   ✅ Matriz 3D cargada exitosamente")
        print(f"   Shape final: {matrix_3d.shape}")
        
        # Información de retorno
        info = {
            'filename': filename,
            'dimensions': (nx, ny, nz),
            'data_type': data_type,
            'is_complex': bool(is_complex),
            'total_elements': total_elements,
            'file_size_mb': os.path.getsize(filename) / (1024*1024)
        }
        
        return matrix_3d, (nx, ny, nz), info

def load_grid_coordinates(base_filename):
    """
    Carga coordenadas del grid desde archivos CSV asociados.
    
    Args:
        base_filename: Archivo base (ej: "output_3D_matrix.bin")
        
    Returns:
        tuple: (grid_x, grid_y, grid_z) o (None, None, None) si no existen
    """
    
    base_path = os.path.splitext(base_filename)[0]
    
    coord_files = {
        'x': f"{base_path}_grid_x.csv",
        'y': f"{base_path}_grid_y.csv", 
        'z': f"{base_path}_grid_z.csv"
    }
    
    coordinates = {}
    
    for coord, filepath in coord_files.items():
        if os.path.exists(filepath):
            try:
                # Leer CSV (asumiendo una sola línea con valores separados por comas)
                with open(filepath, 'r') as f:
                    line = f.readline().strip()
                    values = [float(x.strip()) for x in line.split(',') if x.strip()]
                    coordinates[coord] = np.array(values)
                    print(f"   📍 Coordenadas {coord.upper()}: {len(values)} puntos [{values[0]:.2f}, {values[-1]:.2f}]")
            except Exception as e:
                print(f"   ⚠️  Error leyendo {filepath}: {e}")
                coordinates[coord] = None
        else:
            print(f"   ⚠️  Archivo no encontrado: {filepath}")
            coordinates[coord] = None
    
    # Retornar coordenadas o None si alguna falta
    if all(v is not None for v in coordinates.values()):
        return coordinates['x'], coordinates['y'], coordinates['z']
    else:
        return None, None, None

def analyze_3D_matrix(matrix_3d, grid_x=None, grid_y=None, grid_z=None):
    """
    Analiza estadísticamente la matriz 3D.
    """
    
    print("\n📊 === ANÁLISIS ESTADÍSTICO ===")
    
    # Calcular magnitudes
    magnitudes = np.abs(matrix_3d)
    
    # Estadísticas básicas
    non_zero_mask = magnitudes > 1e-15
    non_zero_count = np.sum(non_zero_mask)
    total_elements = matrix_3d.size
    
    print(f"Elementos totales: {total_elements:,}")
    print(f"Elementos no-cero: {non_zero_count:,} ({100*non_zero_count/total_elements:.2f}%)")
    
    if non_zero_count > 0:
        non_zero_mags = magnitudes[non_zero_mask]
        
        print(f"Magnitud mínima: {np.min(non_zero_mags):.6f}")
        print(f"Magnitud máxima: {np.max(non_zero_mags):.6f}")
        print(f"Magnitud promedio: {np.mean(non_zero_mags):.6f}")
        print(f"Magnitud mediana: {np.median(non_zero_mags):.6f}")
        print(f"Desviación estándar: {np.std(non_zero_mags):.6f}")
        
        # Encontrar posición del máximo
        max_idx = np.unravel_index(np.argmax(magnitudes), magnitudes.shape)
        max_value = matrix_3d[max_idx]
        
        print(f"\n🎯 Máximo global:")
        print(f"   Posición índice: {max_idx}")
        print(f"   Valor complejo: {max_value:.6f}")
        print(f"   Magnitud: {np.abs(max_value):.6f}")
        
        if grid_x is not None and grid_y is not None and grid_z is not None:
            real_pos = (grid_x[max_idx[0]], grid_y[max_idx[1]], grid_z[max_idx[2]])
            print(f"   Posición real: ({real_pos[0]:.3f}, {real_pos[1]:.3f}, {real_pos[2]:.3f}) m")
    
    else:
        print("⚠️  No se encontraron elementos no-cero")
    
    return {
        'non_zero_count': non_zero_count,
        'total_elements': total_elements,
        'magnitudes': magnitudes,
        'non_zero_mask': non_zero_mask
    }

def visualize_3D_matrix(matrix_3d, analysis, grid_x=None, grid_y=None, grid_z=None, output_dir=None):
    """
    Genera visualizaciones de la matriz 3D.
    """
    
    print("\n🎨 === GENERANDO VISUALIZACIONES ===")
    
    magnitudes = analysis['magnitudes']
    non_zero_mask = analysis['non_zero_mask']
    nx, ny, nz = matrix_3d.shape
    
    # Configurar matplotlib
    plt.style.use('default')
    fig = plt.figure(figsize=(20, 12))
    
    # 1. Histograma de magnitudes
    ax1 = plt.subplot(2, 3, 1)
    if analysis['non_zero_count'] > 0:
        non_zero_mags = magnitudes[non_zero_mask]
        plt.hist(non_zero_mags, bins=50, alpha=0.7, color='blue', edgecolor='black')
        plt.xlabel('Magnitud')
        plt.ylabel('Frecuencia')
        plt.title('Distribución de Magnitudes (valores no-cero)')
        plt.yscale('log')
        plt.grid(True, alpha=0.3)
    else:
        plt.text(0.5, 0.5, 'No hay valores\nno-cero', ha='center', va='center', 
                transform=ax1.transAxes, fontsize=14)
        plt.title('Distribución de Magnitudes')
    
    # 2. Proyección XY (suma en Z)
    ax2 = plt.subplot(2, 3, 2)
    xy_projection = np.sum(magnitudes, axis=2)
    if grid_x is not None and grid_y is not None:
        extent = [grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]]
        im2 = plt.imshow(xy_projection.T, origin='lower', aspect='auto', 
                        cmap='hot', extent=extent)
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
    else:
        im2 = plt.imshow(xy_projection.T, origin='lower', aspect='auto', cmap='hot')
        plt.xlabel('Índice X')
        plt.ylabel('Índice Y')
    plt.title('Proyección XY (suma en Z)')
    plt.colorbar(im2, ax=ax2, label='Magnitud acumulada')
    
    # 3. Proyección XZ (suma en Y) 
    ax3 = plt.subplot(2, 3, 3)
    xz_projection = np.sum(magnitudes, axis=1)
    if grid_x is not None and grid_z is not None:
        extent = [grid_x[0], grid_x[-1], grid_z[0], grid_z[-1]]
        im3 = plt.imshow(xz_projection.T, origin='lower', aspect='auto', 
                        cmap='hot', extent=extent)
        plt.xlabel('X (m)')
        plt.ylabel('Z (m)')
    else:
        im3 = plt.imshow(xz_projection.T, origin='lower', aspect='auto', cmap='hot')
        plt.xlabel('Índice X')
        plt.ylabel('Índice Z')
    plt.title('Proyección XZ (suma en Y)')
    plt.colorbar(im3, ax=ax3, label='Magnitud acumulada')
    
    # 4. Proyección YZ (suma en X)
    ax4 = plt.subplot(2, 3, 4)
    yz_projection = np.sum(magnitudes, axis=0)
    if grid_y is not None and grid_z is not None:
        extent = [grid_y[0], grid_y[-1], grid_z[0], grid_z[-1]]
        im4 = plt.imshow(yz_projection.T, origin='lower', aspect='auto', 
                        cmap='hot', extent=extent)
        plt.xlabel('Y (m)')
        plt.ylabel('Z (m)')
    else:
        im4 = plt.imshow(yz_projection.T, origin='lower', aspect='auto', cmap='hot')
        plt.xlabel('Índice Y')
        plt.ylabel('Índice Z')
    plt.title('Proyección YZ (suma en X)')
    plt.colorbar(im4, ax=ax4, label='Magnitud acumulada')
    
    # 5. Corte en Z medio
    ax5 = plt.subplot(2, 3, 5)
    z_mid = nz // 2
    z_slice = magnitudes[:, :, z_mid]
    if grid_x is not None and grid_y is not None:
        extent = [grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]]
        im5 = plt.imshow(z_slice.T, origin='lower', aspect='auto', 
                        cmap='hot', extent=extent)
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        z_val = grid_z[z_mid] if grid_z is not None else z_mid
        plt.title(f'Corte XY en Z={z_val:.2f}')
    else:
        im5 = plt.imshow(z_slice.T, origin='lower', aspect='auto', cmap='hot')
        plt.xlabel('Índice X')
        plt.ylabel('Índice Y')
        plt.title(f'Corte XY en Z={z_mid}')
    plt.colorbar(im5, ax=ax5, label='Magnitud')
    
    # 6. Visualización 3D de puntos no-cero
    ax6 = fig.add_subplot(2, 3, 6, projection='3d')
    if analysis['non_zero_count'] > 0 and analysis['non_zero_count'] < 1000:  # Limitar puntos
        indices = np.where(non_zero_mask)
        
        if grid_x is not None and grid_y is not None and grid_z is not None:
            # Usar coordenadas reales
            x_coords = grid_x[indices[0]]
            y_coords = grid_y[indices[1]]
            z_coords = grid_z[indices[2]]
            ax6.set_xlabel('X (m)')
            ax6.set_ylabel('Y (m)')
            ax6.set_zlabel('Z (m)')
        else:
            # Usar índices
            x_coords = indices[0]
            y_coords = indices[1]
            z_coords = indices[2]
            ax6.set_xlabel('Índice X')
            ax6.set_ylabel('Índice Y')
            ax6.set_zlabel('Índice Z')
        
        colors = magnitudes[non_zero_mask]
        scatter = ax6.scatter(x_coords, y_coords, z_coords, 
                             c=colors, cmap='hot', alpha=0.6, s=20)
        ax6.set_title('Puntos no-cero en 3D')
        plt.colorbar(scatter, ax=ax6, label='Magnitud', shrink=0.5)
    else:
        ax6.text(0.5, 0.5, 0.5, 'Demasiados puntos\npara visualizar', 
                ha='center', va='center', transform=ax6.transAxes)
        ax6.set_title('Visualización 3D (limitada)')
    
    plt.tight_layout()
    
    # Guardar figura si se especifica directorio
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, '3D_matrix_analysis.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"   💾 Visualización guardada: {output_file}")
    
    plt.show()
    print("   ✅ Visualizaciones generadas")

def main():
    """
    Función principal del script.
    """
    
    print("🎯 === VISUALIZADOR DE MATRIZ 3D SAR ===\n")
    
    # Determinar archivo de entrada
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        # Usar archivo por defecto
        default_file = "build/output/proc/output_3D_matrix.bin"
        if os.path.exists(default_file):
            filename = default_file
        else:
            print("❌ Error: Especifique un archivo de matriz 3D")
            print("   Uso: python visualize_3D_matrix.py [archivo.bin]")
            sys.exit(1)
    
    try:
        # Leer matriz 3D
        matrix_3d, dimensions, info = read_3D_complex_matrix(filename)
        
        print(f"\n📋 Información del archivo:")
        for key, value in info.items():
            if key == 'file_size_mb':
                print(f"   {key}: {value:.2f} MB")
            else:
                print(f"   {key}: {value}")
        
        # Cargar coordenadas del grid
        print(f"\n📍 Cargando coordenadas del grid...")
        grid_x, grid_y, grid_z = load_grid_coordinates(filename)
        
        # Analizar matriz
        analysis = analyze_3D_matrix(matrix_3d, grid_x, grid_y, grid_z)
        
        # Generar visualizaciones
        output_dir = os.path.dirname(filename)
        visualize_3D_matrix(matrix_3d, analysis, grid_x, grid_y, grid_z, output_dir)
        
        print(f"\n✅ Análisis completado exitosamente")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()