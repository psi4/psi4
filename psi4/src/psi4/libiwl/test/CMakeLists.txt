# Standalone round-trip exerciser for libiwl. Links the libiwl/libpsio/
# libpsi4util/liboptions/libciomr static archives plus a tiny stubs TU
# (test/stubs.cc) that supplies the runtime globals normally provided by
# core.cc / libdpd, and one no-op leaf (C_DSYEV) that libciomr's unity object
# references but the test never exercises. This keeps the test self-contained
# without pulling the full libmints/Libint closure.
add_executable(iwl_roundtrip_test.x iwl_roundtrip_test.cc stubs.cc)
target_link_libraries(iwl_roundtrip_test.x PRIVATE iwl psio psi4util options ciomr)
target_include_directories(iwl_roundtrip_test.x PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
)
add_test(NAME iwl_roundtrip
         COMMAND iwl_roundtrip_test.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(iwl_roundtrip PROPERTIES LABELS "libiwl;quicktests")
