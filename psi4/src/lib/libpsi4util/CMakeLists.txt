set(sources_list quad.cc 
                 stl_string.cc 
                 memory_manager.cc 
                 exception.cc 
                 combinations.cc 
)
psi4_add_library(psi4util sources_list)
target_link_libraries(psi4util options)
target_include_directories(psi4util PUBLIC options)
