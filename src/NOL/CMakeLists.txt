AUX_SOURCE_DIRECTORY(. DIR_LIB_SRCS)
add_library(${OptimizationLibraryName} STATIC  ${DIR_LIB_SRCS} steepest_descent.h newton.h global.h line_search.h)
target_link_libraries (${OptimizationLibraryName} Eigen3::Eigen)
