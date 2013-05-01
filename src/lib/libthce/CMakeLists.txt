set(SRC laplace.cc lreri.cc   thce.cc    thcew.cc)
add_library(thce ${SRC})
add_dependencies(thce mints)
