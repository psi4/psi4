# Standalone round-trip exerciser for libpsio. Links the libpsio/libpsi4util/
# liboptions/libciomr static archives plus a tiny stubs TU (test/stubs.cc) for
# the runtime globals normally provided by core.cc, so the test is
# self-contained (no python, no libdpd).
add_executable(psio_roundtrip_test.x psio_roundtrip_test.cc stubs.cc)
target_link_libraries(psio_roundtrip_test.x PRIVATE psio psi4util options ciomr)
target_include_directories(psio_roundtrip_test.x PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
)
add_test(NAME psio_roundtrip
         COMMAND psio_roundtrip_test.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(psio_roundtrip PROPERTIES LABELS "libpsio;quicktests")
