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
# globals normally provided by core.cc / libdpd. That curated link resolves
# with GNU/Clang on Linux/macOS but not with Windows lld-link (the unity/PCH
# build pulls in extra symbol references the link line lacks), so restrict it
# to non-Windows rather than drag in the full dependency closure.
if(NOT WIN32)
    add_subdirectory(test)
endif()
