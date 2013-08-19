set(SRC apps.cc cubature.cc direct_screening.cc hamiltonian.cc jk.cc link.cc points.cc solver.cc v.cc wrapper.cc)
add_library(fock ${SRC})
add_dependencies(fock mints)
