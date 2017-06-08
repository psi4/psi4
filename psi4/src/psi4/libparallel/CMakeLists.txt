set(sources_list
                 PsiOutStream.cc 
                 BasesBase.cc 
                 process.cc 
                 ParallelPrinter.cc 
)
psi4_add_module(lib parallel sources_list options)
