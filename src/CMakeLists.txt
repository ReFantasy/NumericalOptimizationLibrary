AUX_SOURCE_DIRECTORY(. DIR_LIB_SRCS)
add_library(opti STATIC  ${DIR_LIB_SRCS} sd.h newton.h global.h line_search.h)
target_link_libraries (opti Eigen3::Eigen)
