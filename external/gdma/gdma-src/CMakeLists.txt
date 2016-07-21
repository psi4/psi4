cmake_minimum_required(VERSION 3.0)
project(gdma LANGUAGES Fortran C)

set(sources_list atom_grids.f90 
                 dma.f90 
                 gdma.f90 
                 input.F90  
                 timing.F90
)
add_library(gdma ${sources_list})

include(FortranCInterface)
FortranCInterface_HEADER(GDMA_MANGLE.h MACRO_NAMESPACE "GDMA_MANGLE_"
                                       SYMBOL_NAMESPACE "GDMA_MANGLE_"
)
install(FILES "${CMAKE_BINARY_DIR}/GDMA_MANGLE.h" 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/GDMA)

install(TARGETS gdma DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

configure_file("gdmaConfig.cmake.in" "gdmaConfig.cmake" @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/gdmaConfig.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX}
)


