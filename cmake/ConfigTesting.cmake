include(ConfigFramework)

# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_CXX_COMPILER_ID}-${BLAS_TYPE}-${CMAKE_BUILD_TYPE}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest own timeout
set(DART_TESTING_TIMEOUT
    "1800"
    CACHE STRING
    "Set timeout in seconds for every single test"
    )

include(CTest)
set(DASHBOARD_DIR ${CMAKE_SOURCE_DIR}/cdash)
add_subdirectory(cdash)
enable_testing()
# This must come last!!
add_subdirectory(tests)
