cmake_minimum_required(VERSION 3.0)
project(Optimization)

set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

include_directories("E:/eigen-3.4.0")


AUX_SOURCE_DIRECTORY(. SRC_LIST)
add_executable(opti ${SRC_LIST} pqf_by_sd.h example.h newton.h global.h)

