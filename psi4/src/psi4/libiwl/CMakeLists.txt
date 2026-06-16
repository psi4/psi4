list(APPEND sources
  buf_close.cc
  buf_fetch.cc
  buf_flush.cc
  buf_init.cc
  buf_put.cc
  buf_wrt.cc
  buf_wrt_mat.cc
  buf_wrt_val.cc
  iwl_reader.cc
  rdone.cc
  wrtone.cc
  )
psi4_add_module(lib iwl sources)

# test/ contains a standalone round-trip exerciser (iwl_roundtrip_test.cc).
# It links against the libiwl/libpsio/libpsi4util/liboptions/libciomr static
# archives plus a tiny stubs TU (test/stubs.cc) that supplies the runtime
# globals normally provided by core.cc / libdpd. That curated link resolves with
# GNU/Clang on Linux/macOS but not with clang-cl/lld-link, which emits inline and
# implicit functions that GCC discards as unused, so the Windows objects carry
# undefined references the link line does not provide: psio.lib wants
# DPDMOSpace::~DPDMOSpace (instantiated by std::vector<DPDMOSpace> in dpd.h), and
# psi4util.lib wants Molecule::Molecule(const Molecule&) (from inline
# Molecule::clone) and ExternalPotential::print (from inline py_print). Not a
# unity or PCH artifact -- Linux is unity-built too, and Azure Windows builds with
# PCH off. Restrict the test to non-Windows rather than drag in the full closure.
if(NOT WIN32)
    add_subdirectory(test)
endif()
