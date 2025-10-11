#include "project_paths.h"
#include <iostream>

// Platform-specific directory creation
#ifdef _WIN32
    #include <direct.h>
    #include <io.h>
    #define mkdir(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <unistd.h>
    #define mkdir(path) mkdir(path, 0755)
#endif

// Helper function to create directory recursively
void createDirectoryRecursive(const std::string& path) {
    std::string current = "";
    for (size_t i = 0; i < path.length(); ++i) {
        if (path[i] == '/' || path[i] == '\\') {
            if (!current.empty()) {
                mkdir(current.c_str());
            }
            current += path[i];
        } else {
            current += path[i];
        }
    }
    mkdir(current.c_str());
}

// Initialize static paths - all centralized under build/
const std::string ProjectPaths::BUILD_DIR = "build";
const std::string ProjectPaths::OUTPUT_DIR = BUILD_DIR + "/output";
const std::string ProjectPaths::BIN_DIR = BUILD_DIR + "/bin";
const std::string ProjectPaths::DATA_DIR = OUTPUT_DIR + "/data";
const std::string ProjectPaths::ASSETS_DIR = "assets";

// Subdirectories under output
const std::string ProjectPaths::PROC_DIR = OUTPUT_DIR + "/proc";
const std::string ProjectPaths::PLOT_DIR = OUTPUT_DIR + "/plot_data";
const std::string ProjectPaths::REPORT_DIR = OUTPUT_DIR + "/reports";

std::string ProjectPaths::getOutputPath(const std::string& filename) {
    return OUTPUT_DIR + "/" + filename;
}

std::string ProjectPaths::getProcPath(const std::string& filename) {
    return PROC_DIR + "/" + filename;
}

std::string ProjectPaths::getPlotPath(const std::string& filename) {
    return PLOT_DIR + "/" + filename;
}

std::string ProjectPaths::getAssetsPath(const std::string& filename) {
    return ASSETS_DIR + "/" + filename;
}

void ProjectPaths::ensureDirectoriesExist() {
    try {
        createDirectoryRecursive(OUTPUT_DIR);
        createDirectoryRecursive(PROC_DIR);
        createDirectoryRecursive(PLOT_DIR);
        createDirectoryRecursive(REPORT_DIR);
        createDirectoryRecursive(DATA_DIR);
        
        std::cout << "[INFO] Created output directories under: " << OUTPUT_DIR << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Failed to create directories: " << e.what() << std::endl;
    }
}