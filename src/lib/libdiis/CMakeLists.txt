set(SRC diisentry.cc diismanager.cc)
add_library(diis ${SRC})
add_dependencies(diis mints)
