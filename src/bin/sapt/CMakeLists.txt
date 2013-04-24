set(SRC wrapper.cc)
add_library(sapt ${SRC})
add_dependencies(sapt mints)
