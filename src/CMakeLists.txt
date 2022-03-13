add_subdirectory(NOL)



# 是否开启测试
option (BUILD_TEST "UnitTest" OFF)
if(BUILD_TEST)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL         ${CMAKE_SOURCE_DIR}/extern/googletest-release-1.11.0.zip
  #GIT_REPOSITORY https://github.com/google/googletest.git
  #GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929
)
FetchContent_MakeAvailable(googletest)
if(WIN32)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
ENDIF()
add_subdirectory(test)
ENDIF()

