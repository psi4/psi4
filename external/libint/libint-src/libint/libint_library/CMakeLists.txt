cmake_minimum_required(VERSION 3.0)
project(libint_library LANGUAGES C CXX)

file(GLOB FULL_SRC "${CMAKE_BINARY_DIR}/*.cc")
file(GLOB FULL_HEAD "${CMAKE_BINARY_DIR}/*.h")
add_library(int STATIC ${FULL_SRC})

install(FILES ${FULL_HEAD}  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/libint)
install(TARGETS int ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
