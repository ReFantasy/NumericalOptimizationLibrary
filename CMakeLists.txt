cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE Release)

find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/lib/eigen/share/eigen3/cmake)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)


AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(main  ${SRC_LIST} example.h )
target_link_libraries (main opti)

# 是否开启测试
option (BUILD_TESTING "Use provided math implementation" OFF)
if (BUILD_TESTING)
add_subdirectory(${CMAKE_SOURCE_DIR}/lib/googletest)
add_subdirectory(test)
ENDIF()
