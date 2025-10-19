cmake_minimum_required(VERSION 3.10)
project(RadarTrajectory)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP for parallel processing
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found. Parallel processing enabled.")
endif()

# Option to use GDAL (set to OFF for simplified version)
option(USE_GDAL "Use GDAL library for real GeoTIFF support" OFF)

# Include directories
include_directories(include)

if(USE_GDAL)
    # Find GDAL package
    find_package(GDAL REQUIRED)
    include_directories(${GDAL_INCLUDE_DIRS})
    
    # Add executable with GDAL support
    add_executable(gs_test_script 
        gs_test_script.cpp
        src/funcao_espiral.cpp
        src/load_dem.cpp
        src/plot_escenario.cpp
        src/crear_grid.cpp
        src/raw_data_processing.cpp
        src/data_export.cpp
        src/fft_utils.cpp
        src/radar_init.cpp
        src/trajectory_generator.cpp
        src/project_paths.cpp
    )
    
    # Add proc_test_script executable with GDAL support
    add_executable(proc_test_script 
        proc_test_script.cpp
        src/funcao_espiral.cpp
        src/radar_init.cpp
        src/trajectory_generator.cpp
        src/binary_reader.cpp
        src/interpolation.cpp
        src/csv_reader.cpp
        src/matrix_reader.cpp
        src/raw_data_processing.cpp
        src/data_export.cpp
        src/fft_utils.cpp
        src/project_paths.cpp
    )
    
    # Link GDAL libraries for proc_test_script
    target_link_libraries(proc_test_script ${GDAL_LIBRARIES})
    target_compile_definitions(proc_test_script PRIVATE USE_GDAL)
    
    message(STATUS "Building with GDAL support")
else()
    # Add executable without GDAL (simplified version)
    add_executable(gs_test_script 
        gs_test_script.cpp
        src/funcao_espiral.cpp
        src/load_dem_simple.cpp
        src/plot_escenario.cpp
        src/crear_grid.cpp
        src/raw_data_processing.cpp
        src/data_export.cpp
        src/fft_utils.cpp
        src/radar_init.cpp
        src/trajectory_generator.cpp
        src/project_paths.cpp
    )
    
    # Add proc_test_script executable without GDAL (simplified version)
    add_executable(proc_test_script 
        proc_test_script.cpp
        src/funcao_espiral.cpp
        src/radar_init.cpp
        src/trajectory_generator.cpp
        src/binary_reader.cpp
        src/interpolation.cpp
        src/csv_reader.cpp
        src/matrix_reader.cpp
        src/raw_data_processing.cpp
        src/data_export.cpp
        src/fft_utils.cpp
        src/project_paths.cpp
    )
    
    message(STATUS "Building simplified version without GDAL")
endif()

# Link math library if needed (for some systems)
if(UNIX)
    target_link_libraries(gs_test_script m)
    target_link_libraries(proc_test_script m)
endif()

# Link OpenMP for parallel processing
if(OpenMP_CXX_FOUND)
    target_link_libraries(gs_test_script OpenMP::OpenMP_CXX)
    target_link_libraries(proc_test_script OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP linked to executables")
endif()

# Set output directory - centralize everything in build/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Ensure consistent output directories for all configurations
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/lib)
endforeach()

# Set output directory
set_target_properties(gs_test_script PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
set_target_properties(proc_test_script PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)