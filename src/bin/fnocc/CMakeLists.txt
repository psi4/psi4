set(SRC blas.cc ccsd.cc density_fitting.cc diis.cc fnocc.cc frozen_natural_orbitals.cc linear.cc local_triples.cc lowmemory_triples.cc mp2.cc opdm.cc quadratic.cc sortintegrals.cc triples.cc)
add_library(fnocc ${SRC})
add_dependencies(fnocc mints)
