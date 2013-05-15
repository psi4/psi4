set(SRC main.cc)
add_library(scf ${SRC})
add_dependencies(scf mints)
