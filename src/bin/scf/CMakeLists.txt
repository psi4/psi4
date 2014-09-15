set(SRC main.cc scf_dummy.cc)
add_library(scf ${SRC})
add_dependencies(scf mints)
