cmake_minimum_required(VERSION 3.16)

project(CppBOLOS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/Eigen)

file(GLOB SOURCES "src/*.cpp")

# Specify the main executable
add_executable(CppBOLOS main.cpp ${SOURCES})


