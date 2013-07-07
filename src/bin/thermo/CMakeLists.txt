set(SRC thermo.cc)
add_library(thermo ${SRC})
add_dependencies(thermo mints)
