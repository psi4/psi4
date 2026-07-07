list(APPEND sources
  mrcc.cc
  )
psi4_add_module(bin mrcc sources)
target_link_libraries(mrcc
  PRIVATE
    pybind11::module
  )
