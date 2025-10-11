@echo off
echo Cleaning up project structure...

REM Remove duplicate folders that get created in wrong places
if exist "bin" rd /s /q "bin"
if exist "output" rd /s /q "output"  
if exist "plot_data" rd /s /q "plot_data"

REM Clean build cache
if exist "CMakeFiles" rd /s /q "CMakeFiles"
if exist "CMakeCache.txt" del "CMakeCache.txt"
if exist "cmake_install.cmake" del "cmake_install.cmake"
if exist "build.ninja" del "build.ninja"
if exist ".ninja_deps" del ".ninja_deps"
if exist ".ninja_log" del ".ninja_log"

echo Creating clean build directory...
if not exist "build" mkdir "build"
cd build

echo Configuring project...
cmake ..

echo Building project...
ninja

echo.
echo =====================================
echo Build completed!
echo Executables: build/bin/
echo Output data: build/output/
echo =====================================

cd ..