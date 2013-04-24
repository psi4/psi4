set(SRC main.cc)
add_library(deriv_wrapper ${SRC})
add_dependencies(deriv_wrapper mints)
