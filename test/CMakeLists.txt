enable_testing()

AUX_SOURCE_DIRECTORY(. TEST_SRCS)

add_executable(
  unit_test
  ${TEST_SRCS}
  test_functions.h)
target_link_libraries(
  unit_test
  gtest_main
  ${OptimizationLibraryName}
)
include(GoogleTest)
gtest_discover_tests(unit_test)