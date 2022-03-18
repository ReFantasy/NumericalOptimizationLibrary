cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

# 库文件输出目录
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
# 可执行文件输出目录
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE Release)

set(OptimizationLibraryName nol)

find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/extern/eigen/share/eigen3/cmake)

#include_directories(${CMAKE_SOURCE_DIR}/src/NOL)
include_directories(${CMAKE_SOURCE_DIR}/src/NOL/include/)
add_subdirectory(src)

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(main  ${SRC_LIST} example.h )
target_link_libraries (main ${OptimizationLibraryName})
