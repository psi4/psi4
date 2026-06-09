list(APPEND sources
  buf_close.cc
  buf_fetch.cc
  buf_flush.cc
  buf_init.cc
  buf_put.cc
  buf_wrt.cc
  buf_wrt_mat.cc
  buf_wrt_val.cc
  rdone.cc
  wrtone.cc
  )
psi4_add_module(lib iwl sources)

# test/ contains a standalone round-trip exerciser (iwl_roundtrip_test.cc)
# that covers the bucket-boundary cases the 1995 README warned about. It is
# not built by default because libpsio's static lib is not standalone --
# PSIO::close references libdpd's `global_dpd_`, and libpsi4util's Options
# system pulls in further hosting-context globals (`psi_file_prefix`,
# `restart_id`) normally provided by the core python bindings TU. Re-enable
# once a proper unit-test harness with stubs / gtest infra lands.
# add_subdirectory(test)
