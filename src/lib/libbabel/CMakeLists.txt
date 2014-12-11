set(libbabel_srcs LibBabel.h
                  LibBabel.cc
                  Psi4OBUnitCell.h
                  Psi4OBUnitCell.cc
)
add_library(babel ${libbabel_srcs})
if(NOT BABEL_FOUND)
  add_dependencies(babel OpenBabel)
endif()
include_directories(${BABEL_INCLUDE_DIR})
target_link_libraries(babel ${BABEL_LIBRARY_DIR}/libopenbabel.so)
