set(SRC mints.cc)
add_library(mints_wrapper ${SRC})
add_dependencies(mints_wrapper mints)
