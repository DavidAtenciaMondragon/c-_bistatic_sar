#!/usr/bin/env python3
"""
Script para visualizar archivos binarios complejos del radar SAR
Archivos soportados: auxData_complex.bin, rawData_complex.bin, rootData_complex.bin

Formato de datos:
- Dimensiones: 256 x 6689 (samples x positions)
- Tipo: double complex (8 bytes real + 8 bytes imaginario por muestra)
- Orden: Real-Imaginario intercalado (R-I-R-I...)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

def read_complex_binary(filepath, rows=None, cols=None):
    """
    Lee archivo binario complejo con el formato específico del radar SAR
    
    Formato del archivo:
    - Header:
      * numDimensions (uint32): Número de dimensiones (2)
      * isComplex (uint32): Flag de complejidad (1=complejo, 0=real)
      * dataTypeLen (uint32): Longitud del string del tipo
      * dataTypeStr (char[16]): Tipo de dato ("double") padded con nulls
      * dimensions (uint64[]): Dimensiones [rows, cols]
    - Data:
      * Valores intercalados R-I-R-I... (double, double, ...)
    
    Args:
        filepath: Ruta al archivo binario
        rows: Número de filas (opcional, se lee del header)
        cols: Número de columnas (opcional, se lee del header)
    
    Returns:
        numpy.ndarray: Matriz compleja de dimensiones (rows, cols)
    """
    try:
        with open(filepath, 'rb') as file:
            print(f"Archivo: {filepath}")
            
            # 1. Leer número de dimensiones (uint32)
            num_dimensions = np.frombuffer(file.read(4), dtype=np.uint32)[0]
            print(f"Número de dimensiones: {num_dimensions}")
            
            # 2. Leer flag de complejidad (uint32)
            is_complex_flag = np.frombuffer(file.read(4), dtype=np.uint32)[0]
            is_complex = (is_complex_flag == 1)
            print(f"Es complejo: {is_complex}")
            
            # 3. Leer longitud del tipo de dato (uint32)
            data_type_len = np.frombuffer(file.read(4), dtype=np.uint32)[0]
            print(f"Longitud del tipo: {data_type_len}")
            
            # 4. Leer tipo de dato (char[16])
            data_type_bytes = file.read(16)
            data_type = data_type_bytes[:data_type_len].decode('utf-8')
            print(f"Tipo de dato: '{data_type}'")
            
            # 5. Leer dimensiones (uint64[])
            dimensions = []
            for i in range(num_dimensions):
                dim = np.frombuffer(file.read(8), dtype=np.uint64)[0]
                dimensions.append(dim)
            
            rows, cols = dimensions[0], dimensions[1]
            print(f"Dimensiones: {rows} x {cols}")
            
            # 6. Leer datos
            total_elements = rows * cols
            values_to_read = total_elements * (2 if is_complex else 1)
            
            # Leer como double (8 bytes cada uno)
            data_bytes = file.read(values_to_read * 8)
            print(f"Bytes de datos leídos: {len(data_bytes)}")
            print(f"Valores esperados: {values_to_read}")
            
            if len(data_bytes) != values_to_read * 8:
                print(f"Warning: Tamaño de datos no coincide. Esperado: {values_to_read * 8}, Actual: {len(data_bytes)}")
                return None
            
            # Convertir a array de double
            data_array = np.frombuffer(data_bytes, dtype=np.float64)
            
            if is_complex:
                # Los datos están intercalados: [R1, I1, R2, I2, ...]
                # Separar partes reales e imaginarias
                real_parts = data_array[0::2]  # Índices pares
                imag_parts = data_array[1::2]  # Índices impares
                
                # Crear array complejo
                complex_data = real_parts + 1j * imag_parts
                
                # Reshape a matriz 2D
                matrix = complex_data.reshape(rows, cols)
            else:
                # Datos reales solamente
                matrix = data_array.reshape(rows, cols)
            
            print(f"Matriz cargada: {matrix.shape}")
            print(f"Rango de valores - Real: [{np.real(matrix).min():.3e}, {np.real(matrix).max():.3e}]")
            if is_complex:
                print(f"Rango de valores - Imag: [{np.imag(matrix).min():.3e}, {np.imag(matrix).max():.3e}]")
            
            return matrix
            
    except Exception as e:
        print(f"Error leyendo archivo {filepath}: {e}")
        import traceback
        traceback.print_exc()
        return None

def plot_complex_data(data, title="Datos Complejos", save_path=None):
    """
    Crea visualizaciones comprehensivas de datos complejos
    
    Args:
        data: Matriz compleja numpy
        title: Título para las gráficas
        save_path: Ruta para guardar las imágenes (opcional)
    """
    if data is None:
        print("No hay datos para graficar")
        return
    
    # Crear figura con subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'{title} - Visualización Completa', fontsize=16, fontweight='bold')
    
    # 1. Magnitud
    im1 = axes[0,0].imshow(np.abs(data), aspect='auto', cmap='viridis', origin='lower')
    axes[0,0].set_title('Magnitud |z|')
    axes[0,0].set_xlabel('Posición')
    axes[0,0].set_ylabel('Muestra')
    plt.colorbar(im1, ax=axes[0,0], label='Magnitud')
    
    # 2. Fase
    im2 = axes[0,1].imshow(np.angle(data), aspect='auto', cmap='hsv', origin='lower')
    axes[0,1].set_title('Fase ∠z (radianes)')
    axes[0,1].set_xlabel('Posición')
    axes[0,1].set_ylabel('Muestra')
    plt.colorbar(im2, ax=axes[0,1], label='Fase (rad)')
    
    # 3. Parte Real
    im3 = axes[0,2].imshow(np.real(data), aspect='auto', cmap='RdBu_r', origin='lower')
    axes[0,2].set_title('Parte Real')
    axes[0,2].set_xlabel('Posición')
    axes[0,2].set_ylabel('Muestra')
    plt.colorbar(im3, ax=axes[0,2], label='Real')
    
    # 4. Parte Imaginaria
    im4 = axes[1,0].imshow(np.imag(data), aspect='auto', cmap='RdBu_r', origin='lower')
    axes[1,0].set_title('Parte Imaginaria')
    axes[1,0].set_xlabel('Posición')
    axes[1,0].set_ylabel('Muestra')
    plt.colorbar(im4, ax=axes[1,0], label='Imaginario')
    
    # 5. Magnitud en dB
    magnitude_db = 20 * np.log10(np.abs(data) + 1e-12)  # Evitar log(0)
    im5 = axes[1,1].imshow(magnitude_db, aspect='auto', cmap='plasma', origin='lower')
    axes[1,1].set_title('Magnitud (dB)')
    axes[1,1].set_xlabel('Posición')
    axes[1,1].set_ylabel('Muestra')
    plt.colorbar(im5, ax=axes[1,1], label='dB')
    
    # 6. Espectro de una posición (ejemplo: columna central)
    center_col = data.shape[1] // 2
    spectrum = data[:, center_col]
    sample_indices = np.arange(len(spectrum))
    
    axes[1,2].plot(sample_indices, np.abs(spectrum), 'b-', alpha=0.7, label='Magnitud')
    axes[1,2].set_title(f'Espectro - Posición {center_col}')
    axes[1,2].set_xlabel('Muestra')
    axes[1,2].set_ylabel('Magnitud')
    axes[1,2].grid(True, alpha=0.3)
    axes[1,2].legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Gráfica guardada en: {save_path}")
    
    plt.show()

def plot_comparison(data_dict, save_path=None):
    """
    Compara múltiples archivos de datos en una sola figura
    
    Args:
        data_dict: Diccionario {nombre: matriz_datos}
        save_path: Ruta para guardar la imagen (opcional)
    """
    n_files = len(data_dict)
    if n_files == 0:
        print("No hay datos para comparar")
        return
    
    fig, axes = plt.subplots(n_files, 2, figsize=(12, 4*n_files))
    if n_files == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle('Comparación de Archivos Binarios', fontsize=16, fontweight='bold')
    
    for i, (name, data) in enumerate(data_dict.items()):
        if data is None:
            continue
            
        # Magnitud
        im1 = axes[i,0].imshow(np.abs(data), aspect='auto', cmap='viridis', origin='lower')
        axes[i,0].set_title(f'{name} - Magnitud')
        axes[i,0].set_xlabel('Posición')
        axes[i,0].set_ylabel('Muestra')
        plt.colorbar(im1, ax=axes[i,0])
        
        # Magnitud en dB
        magnitude_db = 20 * np.log10(np.abs(data) + 1e-12)
        im2 = axes[i,1].imshow(magnitude_db, aspect='auto', cmap='plasma', origin='lower')
        axes[i,1].set_title(f'{name} - Magnitud (dB)')
        axes[i,1].set_xlabel('Posición')
        axes[i,1].set_ylabel('Muestra')
        plt.colorbar(im2, ax=axes[i,1])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Comparación guardada en: {save_path}")
    
    plt.show()

def main():
    """
    Visualiza los 3 archivos binarios específicos en subplots:
    - auxData_complex.bin: valor absoluto
    - rawData_complex.bin: parte real  
    - rootData_complex.bin: valor absoluto
    """
    data_dir = Path('build/output/proc')
    
    print("=== Visualizador de Datos Binarios SAR ===")
    print(f"Directorio de datos: {data_dir}")
    
    # Cargar los 3 archivos específicos
    files_to_load = {
        'auxData': ('auxData_complex.bin', 'Magnitud', lambda x: np.abs(x)),
        'rawData': ('rawData_complex.bin', 'Parte Real', lambda x: np.real(x)),
        'rootData': ('rootData_complex.bin', 'Magnitud', lambda x: np.abs(x))
    }
    
    # Cargar datos
    data_dict = {}
    for name, (filename, title, transform) in files_to_load.items():
        file_path = data_dir / filename
        if file_path.exists():
            print(f"\nCargando {filename}...")
            data = read_complex_binary(file_path)
            if data is not None:
                data_dict[name] = (transform(data), title)
            else:
                print(f"Error cargando {filename}")
        else:
            print(f"Archivo no encontrado: {filename}")
    
    if len(data_dict) == 0:
        print("No se pudieron cargar archivos.")
        return
    
    # Crear visualización
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Radar SAR - Datos Procesados', fontsize=16, fontweight='bold')
    
    # Orden específico para los subplots
    plot_order = ['auxData', 'rawData', 'rootData']
    
    for i, name in enumerate(plot_order):
        if name in data_dict:
            data, title = data_dict[name]
            
            # Crear el plot
            im = axes[i].imshow(data, aspect='auto', cmap='viridis', origin='lower')
            axes[i].set_title(f'{name.replace("Data", " Data")} - {title}')
            axes[i].set_xlabel('Posición')
            axes[i].set_ylabel('Muestra')
            
            # Añadir colorbar
            plt.colorbar(im, ax=axes[i], label=title)
        else:
            # Si no se pudo cargar, mostrar mensaje
            axes[i].text(0.5, 0.5, f'{name}\nNo disponible', 
                        ha='center', va='center', transform=axes[i].transAxes,
                        fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
            axes[i].set_title(f'{name.replace("Data", " Data")}')
    
    plt.tight_layout()
    plt.show()
    
    # Mostrar estadísticas
    print("\n=== Estadísticas ===")
    for name in plot_order:
        if name in data_dict:
            data, title = data_dict[name]
            print(f"{name}: {title}")
            print(f"  Forma: {data.shape}")
            print(f"  Rango: [{data.min():.3e}, {data.max():.3e}]")
            print(f"  Media: {data.mean():.3e}")
            print(f"  Std: {data.std():.3e}")
        else:
            print(f"{name}: No disponible")

if __name__ == "__main__":
    main()