cmake_minimum_required(VERSION 3.0)
project(dkh LANGUAGES Fortran C)

#CMake is stupid and won't let us combine the next two commands...
set(src_files dkh2-dkh4_main.F90)
add_library(dkh ${src_files})

include(FortranCInterface)
FortranCInterface_HEADER(DKH_MANGLE.h MACRO_NAMESPACE "DKH_MANGLE_"
                                      SYMBOL_NAMESPACE "DKH_MANGLE_"
)
install(FILES "${CMAKE_BINARY_DIR}/DKH_MANGLE.h" 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/DKH)

install(TARGETS dkh DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

configure_file("dkhConfig.cmake.in" "dkhConfig.cmake" @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/dkhConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX}
)

