set(SRC fprintf.cc p_fprintf.cc parallel.cc process.cc)
add_library(parallel ${SRC})
add_dependencies(parallel mints)
