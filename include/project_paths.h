#ifndef PROJECT_PATHS_H
#define PROJECT_PATHS_H

#include <string>

// Forward declarations to avoid filesystem include in header
namespace std {
    namespace filesystem {
        class path;
    }
}

// Centralized path configuration
class ProjectPaths {
public:
    // Base directories - always relative to build directory
    static const std::string BUILD_DIR;
    static const std::string OUTPUT_DIR;
    static const std::string BIN_DIR;
    static const std::string DATA_DIR;
    static const std::string ASSETS_DIR;
    
    // Subdirectories
    static const std::string PROC_DIR;
    static const std::string PLOT_DIR;
    static const std::string REPORT_DIR;
    
    // Utility functions
    static std::string getOutputPath(const std::string& filename);
    static std::string getProcPath(const std::string& filename);
    static std::string getPlotPath(const std::string& filename);
    static std::string getAssetsPath(const std::string& filename);
    static void ensureDirectoriesExist();
};

#endif // PROJECT_PATHS_H