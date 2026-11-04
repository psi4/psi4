# C++ unit tests for libpsio.
#
# The executables are built inside the psi4-core ExternalProject, which has its
# own binary directory and does not call enable_testing(). Registration has to
# happen here, in the superbuild, because that is where `enable_testing()` runs
# (cmake/ConfigTesting.cmake) and where CI invokes `ctest`.
ExternalProject_Get_Property(psi4-core BINARY_DIR)

# "smoketests" is what actually gets this run: CI drives scarcely anything through
# ctest because the regression tests are doubled by pytest, and a C++ unit test is
# not eligible for pytest (it is not a `psi4 input.dat`). The smoke and plugin
# labels *are* driven through ctest -- Azure Linux runs `ctest -L "plug|smoke"`
# (devtools/scripts/ci_run_test.py), Azure Windows `--label-regex smoke`.
# "quicktests" is kept alongside it, as 13 other tests do, so the GH macOS
# `ctest -L quick` lane keeps running it too.
add_test(NAME libpsio-roundtrip
         COMMAND ${BINARY_DIR}/src/psi4/libpsio/test/psio_roundtrip_test.x${CMAKE_EXECUTABLE_SUFFIX}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(libpsio-roundtrip PROPERTIES LABELS "psi;quicktests;smoketests;libpsio")
