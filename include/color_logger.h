#ifndef COLOR_LOGGER_H
#define COLOR_LOGGER_H

#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#endif

class ColorLogger {
public:
    // Color codes for Windows console
    enum Color {
        BLACK = 0,
        DARK_BLUE = 1,
        DARK_GREEN = 2,
        DARK_CYAN = 3,
        DARK_RED = 4,
        DARK_MAGENTA = 5,
        DARK_YELLOW = 6,
        GRAY = 7,
        DARK_GRAY = 8,
        BLUE = 9,
        GREEN = 10,
        CYAN = 11,
        RED = 12,
        MAGENTA = 13,
        YELLOW = 14,
        WHITE = 15
    };

    static void setColor(Color color) {
#ifdef _WIN32
        HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        if (hConsole != INVALID_HANDLE_VALUE) {
            SetConsoleTextAttribute(hConsole, color);
        }
#endif
    }

    static void resetColor() {
#ifdef _WIN32
        setColor(GRAY);  // Default color
#endif
    }

    // Colored logging functions
    static void logSuccess(const std::string& message) {
        setColor(GREEN);
        std::cout << "[OK] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logError(const std::string& message) {
        setColor(RED);
        std::cout << "[ERROR] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logWarning(const std::string& message) {
        setColor(YELLOW);
        std::cout << "[WARNING] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logInfo(const std::string& message) {
        setColor(CYAN);
        std::cout << "[INFO] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logTimer(const std::string& message) {
        setColor(MAGENTA);
        std::cout << "[TIME] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logProgress(const std::string& message) {
        setColor(WHITE);
        std::cout << "[PROGRESS] ";
        resetColor();
        std::cout << message << std::endl;
    }

    static void logDebug(const std::string& message) {
        setColor(DARK_GRAY);
        std::cout << "[DEBUG] ";
        resetColor();
        std::cout << message << std::endl;
    }
};

// Convenience macros
#define LOG_SUCCESS(msg) ColorLogger::logSuccess(msg)
#define LOG_ERROR(msg) ColorLogger::logError(msg)
#define LOG_WARNING(msg) ColorLogger::logWarning(msg)
#define LOG_INFO(msg) ColorLogger::logInfo(msg)
#define LOG_TIMER(msg) ColorLogger::logTimer(msg)
#define LOG_PROGRESS(msg) ColorLogger::logProgress(msg)
#define LOG_DEBUG(msg) ColorLogger::logDebug(msg)

#endif // COLOR_LOGGER_H