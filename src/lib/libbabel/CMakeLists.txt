set(libbabel_srcs LibBabel.h
                  LibBabel.cc
                  Psi4OBUnitCell.h
                  Psi4OBUnitCell.cc
)
add_library(babel ${libbabel_srcs})
if(NOT BABEL_FOUND)
  add_dependencies(babel openbabel)
  include_directories(${BABEL_INC_DIR})
endif()
target_link_libraries(babel ${BABEL_LIB_DIR}/libopenbabel.so)
