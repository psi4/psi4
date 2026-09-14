# Standalone round-trip exerciser for libpsio. Links the libpsio/libpsi4util/
# liboptions/libciomr static archives plus a tiny stubs TU (test/stubs.cc) for
# the runtime globals normally provided by core.cc, so the test is
# self-contained (no python, no libdpd).
#
# Only the executable is defined here. See ctest registration in `tests/libpsio/`
# (an add_test() in this scope would not register at all: psi4-core is a separate
# ExternalProject that never calls enable_testing()).
add_executable(psio_roundtrip_test.x psio_roundtrip_test.cc stubs.cc)
target_link_libraries(psio_roundtrip_test.x PRIVATE psio psi4util options ciomr)

# Conformance harness for the libspill-backed libpsio: takes the address of every
# free function and PSIO class method Psi4 declares (materialised, since &C::f is
# never null and the compiler would otherwise fold the test away), so a declared-
# but-undefined entry point surfaces as a link error here rather than in the full
# core link. Then drives a real round trip that also asserts two PSIO instances
# don't share state, and that the store lands where PSIOManager says it does.
add_executable(psio_conformance_test.x psio_conformance_test.cc stubs.cc)
target_link_libraries(psio_conformance_test.x PRIVATE psio psi4util options ciomr)
