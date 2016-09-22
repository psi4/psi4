set(sources_list quad.cc 
                 stl_string.cc 
                 memory_manager.cc 
                 exception.cc 
                 combinations.cc 
)
psi4_add_module(lib psi4util sources_list parallel)
