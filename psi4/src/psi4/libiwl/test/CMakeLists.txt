# Standalone round-trip exerciser for libiwl. Links the libiwl/libpsio/
# libpsi4util/liboptions/libciomr static archives plus a tiny stubs TU
# (test/stubs.cc) that supplies the runtime globals normally provided by
# core.cc / libdpd, and one no-op leaf (C_DSYEV) that libciomr's unity object
# references but the test never exercises. This keeps the test self-contained
# without pulling the full libmints/Libint closure.
#
# Only the executable is defined here. See ctest registration in `tests/libiwl/`
# (an add_test() in this scope would not register at all: psi4-core is a separate
# ExternalProject that never calls enable_testing()).
add_executable(iwl_roundtrip_test.x iwl_roundtrip_test.cc stubs.cc)
target_link_libraries(iwl_roundtrip_test.x PRIVATE iwl psio psi4util options ciomr)
