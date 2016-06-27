set(sources_list wrapper.cc 
                 atom_grids.f90 
                 dma.f90 
                 gdma.f90 
                 input.F90  
                 timing.F90
)
add_library(gdma ${sources_list})
target_link_libraries(gdma mints)


