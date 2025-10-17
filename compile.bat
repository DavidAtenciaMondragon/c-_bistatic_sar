@echo off
echo === Radar Project Compilation Script ===

echo Step 1: Checking current directory...
cd /d "C:\Users\David\Desktop\implementacao_c++\core"
echo Current directory: %CD%

echo.
echo Step 2: Going to build directory...
cd build
echo Build directory: %CD%

echo.
echo Step 3: Running CMake build...
cmake --build . --verbose

echo.
echo Step 4: Checking if executables were created...
if exist "bin\proc_test_script.exe" (
    echo ✓ proc_test_script.exe created successfully!
    dir bin\proc_test_script.exe
) else (
    echo ✗ proc_test_script.exe NOT found!
)

if exist "bin\gs_test_script.exe" (
    echo ✓ gs_test_script.exe created successfully!
    dir bin\gs_test_script.exe
) else (
    echo ✗ gs_test_script.exe NOT found!
)

echo.
echo === Compilation completed ===
pause