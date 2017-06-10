set(sources_list stl_string.cc 
                 PsiOutStream.cc
                 process.cc
                 memory_manager.cc 
                 exception.cc 
                 combinations.cc 
)
psi4_add_module(lib psi4util sources_list)
