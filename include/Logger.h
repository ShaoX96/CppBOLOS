// Logger.h
#pragma once

#include <iostream>
namespace CppBOLOS {

enum LogLevel {
    LOG_NONE = 0,
    LOG_WARNING,
    LOG_INFO,
    LOG_DEBUG
};

// Global logging level variable
extern LogLevel currentLogLevel;

#define LOG_DEBUG(msg) if (CppBOLOS::currentLogLevel >= CppBOLOS::LOG_DEBUG) { std::cout << "[BOLOS-DEBUG] " << msg << std::endl; }
#define LOG_INFO(msg) if (CppBOLOS::currentLogLevel >= CppBOLOS::LOG_INFO) { std::cout << "[BOLOS-INFO] " << msg << std::endl; }
#define LOG_WARNING(msg) if (CppBOLOS::currentLogLevel >= CppBOLOS::LOG_WARNING) { std::cerr << "[BOLOS-WARNING] " << msg << std::endl; }

}
