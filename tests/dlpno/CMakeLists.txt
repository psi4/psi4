# C++ unit tests for dlpno.
#
# The executables are built inside the psi4-core ExternalProject, which has its
# own binary directory and does not call enable_testing(). Registration has to
# happen here, in the superbuild, because that is where `enable_testing()` runs
# (cmake/ConfigTesting.cmake) and where CI invokes `ctest`.
ExternalProject_Get_Property(psi4-core BINARY_DIR)

# Labels as for tests/libpsio/: "smoketests" gets this run by Azure Linux
# (`ctest -L "plug|smoke"`) and Azure Windows (`--label-regex smoke`), and
# "quicktests" by the GH macOS `ctest -L quick` lane. A C++ unit test is not
# eligible for pytest.
add_test(NAME dlpno-invert-map
         COMMAND ${BINARY_DIR}/src/psi4/dlpno/test/dlpno_invert_map_test.x${CMAKE_EXECUTABLE_SUFFIX}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(dlpno-invert-map PROPERTIES LABELS "psi;quicktests;smoketests;dlpno")
