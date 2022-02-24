cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_BUILD_TYPE Release)

find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/eigen/share/eigen3/cmake)

#find_package(Ceres REQUIRED PATHS "E:/ceres-solver-2.0.0/install/CMake")

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)

AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(test ${SRC_LIST} example.h )
target_link_libraries (test opti)#Ceres::ceres
