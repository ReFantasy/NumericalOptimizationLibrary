cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

#include_directories("E:/eigen-3.4.0")
#/usr/local/share/eigen3/cmake
find_package (Eigen3 REQUIRED NO_MODULE HINTS "E:/eigen-3.4.0/install/share/eigen3/cmake")

#find_package(Ceres REQUIRED PATHS "E:/ceres-solver-2.0.0/install/CMake")


AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(opti ${SRC_LIST} sd.h example.h newton.h global.h line_search.h)
target_link_libraries (opti Eigen3::Eigen)#Ceres::ceres
