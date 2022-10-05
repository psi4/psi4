list(APPEND sources
  diismanager.cc
  )
psi4_add_module(lib diis sources)

target_link_libraries(diis PUBLIC pybind11::headers)
