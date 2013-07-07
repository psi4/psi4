set(SRC lmp2.cc)
add_library(lmp2 ${SRC})
add_dependencies(lmp2 mints trans dpd)
