set(SRC cuhf.cc frac.cc hf.cc ks.cc mom.cc rhf.cc rohf.cc sad.cc uhf.cc)
add_library(scf_solver ${SRC})
add_dependencies(scf_solver mints)
