cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

set(OptimizationLibraryName OptimizationLibrary)

find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/lib/eigen/share/eigen3/cmake)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(src)

# 是否开启测试
option (BUILD_TEST "UnitTest" OFF)
if(BUILD_TEST)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL         ${CMAKE_SOURCE_DIR}/lib/googletest-release-1.11.0.zip
  #GIT_REPOSITORY https://github.com/google/googletest.git
  #GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929
)
FetchContent_MakeAvailable(googletest)
if(WIN32)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
ENDIF()
add_subdirectory(test)
ENDIF()

AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(main  ${SRC_LIST} example.h )
target_link_libraries (main ${OptimizationLibraryName})