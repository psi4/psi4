# Standalone unit test for invert_map (dlpno/invert_map.h). invert_map is
# header-only, so the test links no psi4 library and needs no stubs TU. The one
# extra source, libpsi4util/exception.cc, provides PsiException, which invert_map
# throws. Compiling it in directly avoids linking the psi4util archive, whose
# unity object would also pull in the core.cc runtime globals.
#
# Only the executable is defined here. See ctest registration in `tests/dlpno/`
# (an add_test() in this scope would not register at all: psi4-core is a separate
# ExternalProject that never calls enable_testing()).
add_executable(dlpno_invert_map_test.x
    invert_map_test.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../../libpsi4util/exception.cc
    )
