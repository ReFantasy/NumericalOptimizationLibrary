cmake_minimum_required(VERSION 3.0)
project(Optimization)

ADD_DEFINITIONS(-DUNICODE) 
ADD_DEFINITIONS(-D_UNICODE) 

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE Release)

#Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

find_package (Eigen3 REQUIRED NO_MODULE HINTS ${CMAKE_SOURCE_DIR}/lib/eigen/share/eigen3/cmake)

include_directories(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(test)
add_subdirectory(src)


AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(main  ${SRC_LIST} example.h )
target_link_libraries (main opti)
